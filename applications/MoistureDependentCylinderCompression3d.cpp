#include "nuto/math/EigenSparseSolve.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"

#include "integrands/GradientDamageAdditive.h"
#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportBoundary.h"
#include "integrands/MoistureTransportCoefficients.h"
#include "structures/MultiPhysicsStructureNew.h"
#include "integrands/Shrinkage.h"

#include <fstream>
#include <iostream>

using namespace NuTo;
using namespace NuTo::Constraint;
using namespace NuTo::Integrands;
using namespace NuTo::Laws;
using namespace NuTo::Visualize;

class OutputFile
{
public:
    OutputFile(std::string fileName)
    {
        mFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
        if (!mFile.is_open())
            throw Exception(__PRETTY_FUNCTION__, "Can't open output file \"" + fileName + "\" !");
    }

    ~OutputFile()
    {
        mFile.close();
    }

    template <typename T>
    OutputFile& Write(const T& data)
    {
        mFile << data;
        mFile.flush();
        return *this;
    }

    template <typename T, typename... TArgs>
    OutputFile& Write(const T& data, const TArgs&... args)
    {
        mFile << data;
        Write(args...);
        return *this;
    }

private:
    std::ofstream mFile;
};

struct CompressionTestSetup
{
    double c = 25.;
    double gf = 0.0045;
    double initialRH = 1.0;
    double equilibriumWV = 0.18;
    double envRelativeHumdity = 1.0;
    double t_dry = 0.;
    double shrinkageRatio = 1.0;
    std::string resultFolder{"MDCCResults"};
};

double TimeDependentDisplacement(double t)
{
    return -0.1 * t;
}

void SaveStressStrain(Group<CellInterface>& groupVolumeCellsTotal, GradientDamage<3, Shrinkage<3>>& integrandVolume,
                      DofType dofDisplacements, double t, const CompressionTestSetup& setup, bool isDrying)
{
    static OutputFile file(setup.resultFolder + "/CylinderStressStrain_c" + std::to_string(setup.c) + "_gf" +
                           std::to_string(setup.gf) + "_rh" + std::to_string(setup.envRelativeHumdity) + "_tdry" +
                           std::to_string(setup.t_dry) + ".dat");

    static double axialStressFromDrying;

    double V{0};
    Eigen::VectorXd Stress{Eigen::VectorXd::Zero(6)};
    for (CellInterface& cell : groupVolumeCellsTotal)
    {
        V += cell.Integrate([&](const CellIpData& cipd) { return 1.0; });
        auto cellStressIntegral = cell.Integrate([&](const CellIpData& cipd) {
            DofVector<double> d;
            //            d[dofDisplacements] = integrandVolume.mLinearElasticDamage.Stress(
            //                    cipd.Apply(dofDisplacements, Nabla::Strain()),
            //                    integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)));
            d[dofDisplacements] = integrandVolume.StressMechanics(cipd);
            return d;
        });
        Stress += cellStressIntegral[dofDisplacements];
    }
    double AxialStress{Stress[2] / V};
    double AxialStrain{TimeDependentDisplacement(t) / 300};

    if (isDrying)
    {
        AxialStrain = 0.;
        axialStressFromDrying = AxialStress;
    }

    file.Write(AxialStrain, ", ", AxialStress - axialStressFromDrying, ", ", AxialStress, "\n");

    std::cout << "Stress (total): " << AxialStress << std::endl;
    std::cout << "Stress (external): " << AxialStress - axialStressFromDrying << std::endl;
    std::cout << "Strain (external): " << AxialStrain << std::endl;
}


void UpdateSetup(CompressionTestSetup& setup, std::string tag, std::string value)
{
    if (tag.compare("C") == 0)
    {
        setup.c = std::stod(value);
        return;
    }
    if (tag.compare("GF") == 0)
    {
        setup.gf = std::stod(value);
        return;
    }
    if (tag.compare("RES") == 0)
    {
        setup.resultFolder = value;
        return;
    }
    if (tag.compare("RHENV") == 0)
    {
        setup.envRelativeHumdity = std::stod(value);
        return;
    }
    if (tag.compare("RHINI") == 0)
    {
        setup.initialRH = std::stod(value);
        return;
    }
    if (tag.compare("TDRY") == 0)
    {
        setup.t_dry = std::stod(value);
        return;
    }
    if (tag.compare("SSC") == 0)
    {
        setup.shrinkageRatio = std::stod(value);
        return;
    }
    if (tag.compare("WVEQ") == 0)
    {
        setup.equilibriumWV = std::stod(value);
        return;
    }
    throw Exception(__PRETTY_FUNCTION__, "Unknown tag: \"" + tag + "\". Use -H for help.");
}


void PrintHelp()
{
    std::cout << "Structure for input values: tag=val" << std::endl;
    std::cout << "Possible tags:" << std::endl;
    std::cout << "  c     : damage law c parameter" << std::endl;
    std::cout << "  gf    : damage law gf parameter" << std::endl;
    std::cout << "  res   : result folder" << std::endl;
    std::cout << "  rhenv : environmental relative humidity [0.0, 1.0]" << std::endl;
    std::cout << "  rhini : initial relative humidity [0.0, 1.0]" << std::endl;
    std::cout << "  ssc   : Stress based shrinkage contribution [0.0, 1.0]" << std::endl;
    std::cout << "  tdry  : drying time [0.0, inf]" << std::endl;
    std::cout << "  wveq  : equilibrium water volume fraction at 100% relative humidity [0.0, inf]" << std::endl;
    exit(EXIT_SUCCESS);
}

std::string FormatString(std::string argument)
{
    std::transform(argument.begin(), argument.end(), argument.begin(), ::toupper);
    return argument;
}


std::pair<std::string, std::string> GetTagAndValue(const std::string& argument)
{
    size_t pos = argument.find_first_of('=', 0);
    if (pos == std::string::npos)
        throw Exception(__PRETTY_FUNCTION__,
                        "Argument does not provide an \"=\" operator. Supported argument structure: tag=val");
    if (pos == 0)
        throw Exception(__PRETTY_FUNCTION__, "No argument tag provided. Supported argument structure: tag=val");
    if (pos == argument.size() - 1)
        throw Exception(__PRETTY_FUNCTION__, "No argument value provided. Supported argument structure: tag=val");

    return {FormatString(argument).substr(0, pos), argument.substr(pos + 1)};
}

void PrintSetup(const CompressionTestSetup& setup)
{
    std::cout << "-------------" << std::endl;
    std::cout << "Program setup" << std::endl;
    std::cout << "-------------" << std::endl;
    std::cout << "Result folder : " << setup.resultFolder << std::endl;
    std::cout << "Damage law c parameter : " << setup.c << std::endl;
    std::cout << "Damage law gf parameter : " << setup.gf << std::endl;
    std::cout << "Environmental relative humidity: " << setup.envRelativeHumdity << std::endl;
    std::cout << "Initial relative humidity: " << setup.initialRH << std::endl;
    std::cout << "Drying time: " << setup.t_dry << std::endl;
    std::cout << "Shrinkage ratio (stress based contribution): " << setup.shrinkageRatio << std::endl;
    std::cout << "Equilibrium water volume fraction at 100% relative humidity: " << setup.equilibriumWV << std::endl;
    std::cout << std::endl << std::endl;
}

CompressionTestSetup SetupTestParameters(int argc, char* argv[])
{
    CompressionTestSetup setup;

    for (int i = 1; i < argc; ++i)
    {
        std::string argument{argv[i]};

        if (FormatString(argument).compare("HELP") == 0 || FormatString(argument).compare("-HELP") == 0 ||
            FormatString(argument).compare("H") == 0 || FormatString(argument).compare("-H") == 0)
            PrintHelp();


        auto tagValuePair = GetTagAndValue(argument);

        UpdateSetup(setup, tagValuePair.first, tagValuePair.second);
    }

    return setup;
}


int main(int argc, char* argv[])
{

    CompressionTestSetup setup = SetupTestParameters(argc, argv);
    PrintSetup(setup);


    std::cout << "Load mesh..." << std::endl;
    auto meshGmsh = MeshGmsh{"Cylinder.msh"};
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};

    // Get element groups
    const auto& groupVolumeElementsTotal = meshGmsh.GetPhysicalGroup("Volume");
    const auto& groupSurfaceElementsUpperFace = meshGmsh.GetPhysicalGroup("UpperFace");
    const auto& groupSurfaceElementsSideFace = meshGmsh.GetPhysicalGroup("SideFace");
    const auto& groupSurfaceElementsLowerFace = meshGmsh.GetPhysicalGroup("LowerFace");
    const auto groupSurfaceElementsTotal =
            Unite(groupSurfaceElementsUpperFace, groupSurfaceElementsSideFace, groupSurfaceElementsLowerFace);

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("displacements", 3);
    ScalarDofType dofNonLocal{"non local equivalent strain"};
    const auto& dofWaterVolumeFraction = MPS.AddDofType("water volume fraction", 1);
    const auto& dofRelativeHumidity = MPS.AddDofType("relative humidity", 1);

    // Create interpolation
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));
    const auto& integrationTriangle5 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofNonLocal, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupSurfaceElementsTotal, interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupSurfaceElementsTotal, interpolationTriangleQuadratic);

    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);


    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsTotal = MPS.CreateCells(groupVolumeElementsTotal, integrationTetrahedron3);
    auto groupSurfaceCellsTotal = MPS.CreateCells(groupSurfaceElementsTotal, integrationTriangle5);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;


    //    auto diffusionWV = MTCoefficientLinear{0., -10., 0., 0., 11.};
    auto diffusionWV = MTCoefficientQuadratic{{{0., 0.}}, {{-10., 0.}}, {{0., 0.}}, {{0., 0.}}, 10.1};
    auto diffusionRH = MTCoefficientConstant{0.};
    auto massExchange = MTCoefficientConstant{10.};
    auto sorptionIso = MTCoefficientLinear{0.0, setup.equilibriumWV, 0., 0., 0.};

    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, diffusionWV, diffusionRH, massExchange,
                              sorptionIso, 1., 0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});


    auto& integrandMoistureTransportBoundary = MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundary, double>{
            MoistureTransportBoundary{integrandMoistureTransport.GetIntegrand(), dofWaterVolumeFraction,
                                      dofRelativeHumidity, 8., 0., setup.envRelativeHumdity},
            &MoistureTransportBoundary::Gradient, &MoistureTransportBoundary::Stiffness});


    Shrinkage<3> lawShrinkage(dofDisplacements, dofWaterVolumeFraction, dofRelativeHumidity, setup.shrinkageRatio,
                              0.000001);

    auto material = Material::DefaultConcrete();
    material.c = setup.c;
    material.gf = setup.gf;
    GradientDamage<3, Shrinkage<3>> integrandVolume{dofDisplacements, dofNonLocal, material, lawShrinkage};
    int nIp = integrationTetrahedron3.GetNumIntegrationPoints();
    integrandVolume.mKappas = Eigen::MatrixXd::Zero(groupVolumeCellsTotal.Size(), nIp);


    // Set nodal values ---------------------------------------------------------------------------

    std::cout << "Set nodal values..." << std::endl;

    MPS.SetNodalValues(dofRelativeHumidity, {setup.initialRH});
    MPS.SetNodalValues(dofWaterVolumeFraction, {sorptionIso.value(0., setup.initialRH, 0., 0.)});
    std::cout << "Nodal values set..." << std::endl;

    // Define constraints -------------------------------------------------------------------------

    Constraints& constraintsDrying = MPS.GetConstraints();


    auto& nodeCenterBottom = mesh.NodeAtCoordinate(Eigen::Vector3d{0., 0., 0.}, dofDisplacements);
    auto& nodeCenterTop = mesh.NodeAtCoordinate(Eigen::Vector3d{0., 0., 300.}, dofDisplacements);
    auto& nodeLeftBottom = mesh.NodeAtCoordinate(Eigen::Vector3d{-50., 0., 0.}, dofDisplacements);
    constraintsDrying.Add(dofDisplacements,
                          {Direction(nodeCenterBottom, Eigen::Vector3d{1., 0., 0.}),
                           Direction(nodeCenterBottom, Eigen::Vector3d{0., 1., 0.}),
                           Direction(nodeCenterBottom, Eigen::Vector3d{0., 0., 1.}),
                           Direction(nodeLeftBottom, Eigen::Vector3d{0., 1., 0.}),
                           Direction(nodeCenterTop, Eigen::Vector3d{0., 1., 0.}),
                           Direction(nodeCenterTop, Eigen::Vector3d{1., 0., 0.})});

    MPS.RenumberDofs();

    // Setup time integration ---------------------------------------------------------------------

    std::cout << "Setup time integration scheme..." << std::endl;

    double t = 0.;
    double dt_max = 0.1;
    double dt_min = dt_max * 1.e-6;
    double dt = dt_max;


    // Newmark

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    // Quasi static
    TimeDependentProblem TDP{&mesh};
    TDP.AddGradientFunction(groupVolumeCellsTotal, [&](const CellIpData& cipd, double t, double dt) {
        return integrandVolume.Gradient(cipd);
    });
    TDP.AddHessian0Function(groupVolumeCellsTotal, [&](const CellIpData& cipd, double t, double dt) {
        return integrandVolume.Hessian0(cipd);
    });
    TDP.AddUpdateFunction(groupVolumeCellsTotal,
                          [&](const CellIpData& cipd, double t, double dt) { return integrandVolume.Update(cipd); });


    // QuasistaticSolver QSS{TDP, dofDisplacements};
    QuasistaticSolver QSS{TDP, {dofDisplacements, dofNonLocal}};
    QSS.SetGlobalTime(t);
    QSS.SetConstraints(constraintsDrying);
    QSS.mTolerance = 1e-4;


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    PostProcess pp(setup.resultFolder);
    pp.DefineVisualizer("Volume", groupVolumeCellsTotal, AverageHandler());
    pp.Add("Volume", dofDisplacements);
    pp.Add("Volume", dofNonLocal);
    pp.Add("Volume", dofWaterVolumeFraction);
    pp.Add("Volume", dofRelativeHumidity);
    pp.Add("Volume",
           [&](const CellIpData& cipd) { return integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)); },
           "Damage");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.Stress(cipd); }, "Stress");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.StressMechanics(cipd); }, "StressMechanics");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.StressAdditive(cipd); }, "StressShrinkage");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.Strain(cipd); }, "Strain");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.StrainMechanics(cipd); }, "StrainMechanics");
    pp.Add("Volume", [&](const CellIpData& cipd) { return integrandVolume.StrainAdditive(cipd); }, "StrainShrinkage");


    // Solve --------------------------------------------------------------------------------------


    if (setup.initialRH < 1.0 && setup.t_dry <= 0.)
    {
        std::cout << "Initial relative humidity is smaller than 1. Calculating mechanical equilibirum" << std::endl;
        QSS.DoStep(0, "MumpsLU");
    }

    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);

    // Drying Process ----------------------------------------------------

    SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t, setup, true);

    double delta_plot = setup.t_dry / 10;
    double next_plot = delta_plot;
    int num_converges = 0;

    std::cout << "Extract nodal values..." << std::endl;

    Eigen::VectorXd d_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 0);
    Eigen::VectorXd v_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 1);
    Eigen::VectorXd a_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 2);
    Eigen::VectorXd delta_d_MT = Eigen::VectorXd::Zero(d_MT.rows());

    std::cout << "Start time integration scheme..." << std::endl;
    EigenSparseSolver solver("MumpsLU");
    while (t < setup.t_dry)
    {
        t += dt;
        std::cout << std::endl
                  << "Next timestep" << std::endl
                  << "-------------" << std::endl
                  << "t  : " << t << std::endl
                  << "dt : " << dt << std::endl;
        // Moisture Transport - Newmark/Crank-Nicholsen

        Eigen::VectorXd d_it = d_MT + delta_d_MT;
        Eigen::VectorXd v_it =
                gamma / (dt * beta) * delta_d_MT + (1 - gamma / beta) * v_MT + dt * (1 - gamma / (2 * beta)) * a_MT;
        Eigen::VectorXd a_it = 1 / (dt * dt * beta) * delta_d_MT - 1 / (dt * beta) * v_MT - (1 / (2 * beta) - 1) * a_MT;

        MPS.MergeDofs(d_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
        MPS.MergeDofs(v_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
        MPS.MergeDofs(a_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);


        Eigen::VectorXd grad = MPS.AssembleResidual({dofWaterVolumeFraction, dofRelativeHumidity},
                                                    {{&groupVolumeCellsTotal, &integrandMoistureTransport},
                                                     {&groupSurfaceCellsTotal, &integrandMoistureTransportBoundary}});

        int iteration = 0;
        while (grad.lpNorm<Eigen::Infinity>() > 10e-13)
        {

            std::cout << "Assemble moisture problem..." << std::endl;
            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K =
                    MPS.AssembleStiffness({dofWaterVolumeFraction, dofRelativeHumidity},
                                          {{&groupVolumeCellsTotal, &integrandMoistureTransport},
                                           {&groupSurfaceCellsTotal, &integrandMoistureTransportBoundary}});

            Eigen::SparseMatrix<double> D =
                    MPS.AssembleDamping({dofWaterVolumeFraction, dofRelativeHumidity},
                                        {{&groupVolumeCellsTotal, &integrandMoistureTransport}});
            Eigen::SparseMatrix<double> H = gamma / (dt * beta) * D + K;

            std::cout << "Solve moisture problem..." << std::endl;
            delta_d_MT = solver.Solve(H, -grad);

            d_it += delta_d_MT;
            v_it += gamma / (dt * beta) * delta_d_MT;
            a_it += 1 / (dt * dt * beta) * delta_d_MT;

            MPS.MergeDofs(d_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
            MPS.MergeDofs(v_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
            MPS.MergeDofs(a_it, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);

            grad = MPS.AssembleResidual({dofWaterVolumeFraction, dofRelativeHumidity},
                                        {{&groupVolumeCellsTotal, &integrandMoistureTransport},
                                         {&groupSurfaceCellsTotal, &integrandMoistureTransportBoundary}});

            std::cout << "Time: " << t << std::endl
                      << "Iteration: " << iteration << std::endl
                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
                      << std::endl;
        }


        std::cout << "Assemble and solve mechanics problem..." << std::endl;
        try
        {
            QSS.DoStep(t, "MumpsLU");
        }
        catch (NewtonRaphson::NoConvergence e)
        {
            std::cout << "No convergence... reducing timestep" << std::endl;
            t -= dt;
            dt *= 0.5;
            num_converges = 0;
            if (dt < dt_min)
                throw e;

            MPS.MergeDofs(d_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
            MPS.MergeDofs(v_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
            MPS.MergeDofs(a_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);
            continue;
        }

        if (std::abs(1. - t) <= 1.e-6)
            dt_max = dt = 1.;
        if (std::abs(2. - t) <= 1.e-6)
            dt_max = dt = 2.;
        if (std::abs(10. - t) <= 1.e-6)
            dt_max = dt = 5.;


        ++num_converges;
        if (dt < dt_max && num_converges > 3)
        {
            std::cout << "Increasing timestep" << std::endl;
            dt *= 2.0;
            if (dt > dt_max)
                dt = dt_max;
            if (t + dt > setup.t_dry)
                dt = setup.t_dry - t;
            num_converges = 0;
        }

        if (t >= next_plot)
        {
            std::cout << std::endl << "Visualize..." << std::endl << std::endl;
            pp.Plot(t, false);
            next_plot += delta_plot;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;
        SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t, setup, true);
    }


    // Compression Test ----------------------------------------------


    t = 0;
    double t_compr = 40;
    dt_max = t_compr / 40;
    next_plot = dt_max;
    delta_plot = dt_max / 10;
    dt = dt_max;

    double dt_max_current = dt_max;
    bool timestepReduced = false;

    double t_adjustBifurcation = dt_max * 3;
    double adjustBifurcation = TimeDependentDisplacement(t_adjustBifurcation) * 0.2;


    Constraints constraintsCompression;

    std::cout << "Apply compression test constraints" << std::endl << std::endl;

    auto groupNodesBottom = mesh.NodesAtAxis(eDirection::Z, dofDisplacements);
    auto groupNodesTop = mesh.NodesAtAxis(eDirection::Z, dofDisplacements, 300);

    std::set<NodeSimple*> handledNodes;
    for (ElementCollectionFem& ec : groupVolumeElementsTotal)
    {
        ElementFem& coordE = ec.CoordinateElement();
        ElementFem& dofE = ec.DofElement(dofDisplacements);
        for (int i = 0; i < dofE.GetNumNodes(); ++i)
        {
            NodeSimple& dofNode = dofE.GetNode(i);
            if (handledNodes.find(&dofNode) != handledNodes.end())
                continue;
            handledNodes.insert(&dofNode);
            Eigen::VectorXd coords = Interpolate(coordE, dofE.Interpolation().GetLocalCoords(i));
            if (coords[2] < 300 - 1e-4 && coords[2] > 0 + 1e-4)
                continue;
            Eigen::VectorXd disp = dofNode.GetValues(0);
            constraintsCompression.Add(dofDisplacements, Direction(dofNode, Eigen::Vector3d{1, 0., 0.}, disp[0]));
            constraintsCompression.Add(dofDisplacements, Direction(dofNode, Eigen::Vector3d{0., 1, 0.}, disp[1]));
            constraintsCompression.Add(dofDisplacements, Direction(dofNode, Eigen::Vector3d{0., 0., 1}, [=](double t) {

                                           double direction = 1.0;
                                           if (coords[2] < 100)
                                               direction = -1;
                                           double dispConst = disp[2] + TimeDependentDisplacement(t) * 0.5 * direction;
                                           double normalizedDistortion =
                                                   ((coords[0] * direction + 50.) / 100.
                                                    // + (std::sin(coords[1] * direction / 50. * 3.14 / 2.) + 1.) / 2.
                                                    ) /
                                                   2.0 * direction;
                                           if (t < t_adjustBifurcation)
                                               dispConst += normalizedDistortion * t / t_adjustBifurcation *
                                                            adjustBifurcation;
                                           else
                                               dispConst += normalizedDistortion * adjustBifurcation;

                                           return dispConst;
                                       }));
        }
    }


    MPS.RenumberDofs();

    QSS.SetGlobalTime(t);
    QSS.SetConstraints(constraintsCompression);


    num_converges = 0;


    std::cout << "Start compression test" << std::endl;
    std::cout << "----------------------" << std::endl << std::endl;

    std::cout << "Start time integration scheme..." << std::endl;
    while (t < t_compr)
    {
        t += dt;
        std::cout << std::endl
                  << "Next timestep" << std::endl
                  << "-------------" << std::endl
                  << "t  : " << t + setup.t_dry << std::endl
                  << "dt : " << dt << std::endl;


        std::cout << "Assemble and solve mechanics problem..." << std::endl;
        try
        {
            QSS.DoStep(t, "MumpsLU");
        }
        catch (NewtonRaphson::NoConvergence e)
        {
            std::cout << "No convergence... reducing timestep" << std::endl;
            t -= dt;
            dt *= 0.125;
            num_converges = 0;
            if (dt < dt_min * 1e-6)
                throw e;

            MPS.MergeDofs(d_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
            MPS.MergeDofs(v_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
            MPS.MergeDofs(a_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);
            continue;
        }

        ++num_converges;
        if (t + dt > 0.0125 * 300 && timestepReduced == false)
        {
            dt_max_current /= 10;
            dt = dt_max_current;
            num_converges = 0;
            timestepReduced = true;
        }
        if (dt < dt_max_current && num_converges > 3)
        {
            std::cout << "Increasing timestep" << std::endl;
            dt *= 2.;
            if (dt > dt_max_current)
                dt = dt_max_current;
            num_converges = 0;
        }

        if (t >= next_plot)
        {
            std::cout << std::endl << "Visualize..." << std::endl << std::endl;
            pp.Plot(t + setup.t_dry, false);
            next_plot = (std::round(t / delta_plot) + 1) * delta_plot;
        }
        SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t, setup, false);
    }
}
