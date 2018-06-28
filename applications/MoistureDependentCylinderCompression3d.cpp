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
    double envRelativeHumdity = 0.4;
    double t_dry = 400.;
    double shrinkageRatio = 0.5;
    std::string resultFolder{"MDCCResults"};
};

double TimeDependentDisplacement(double t)
{
    return -0.1 * t;
}

void SaveStressStrain(Group<CellInterface>& groupVolumeCellsTotal, GradientDamage<3, Shrinkage<3>>& integrandVolume,
                      DofType dofDisplacements, double t, const CompressionTestSetup& setup)
{
    static OutputFile file(setup.resultFolder + "/CylinderStressStrain.dat");
    double V{0};
    Eigen::VectorXd Stress{Eigen::VectorXd::Zero(6)};
    for (CellInterface& cell : groupVolumeCellsTotal)
    {
        V += cell.Integrate([&](const CellIpData& cipd) { return 1.0; });
        auto test = cell.Integrate([&](const CellIpData& cipd) {
            DofVector<double> d;
            d[dofDisplacements] = integrandVolume.mLinearElasticDamage.Stress(
                    cipd.Apply(dofDisplacements, Nabla::Strain()),
                    integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)));
            return d;
        });
        Stress += test[dofDisplacements];
    }
    double AxialStress{Stress[2] / V};
    double AxialStrain{TimeDependentDisplacement(t) / 300};

    file.Write(AxialStress, ", ", AxialStrain, "\n");

    std::cout << "Stress: " << AxialStress << std::endl;
    std::cout << "Strain: " << AxialStrain << std::endl;
}


void UpdateSetup(CompressionTestSetup& setup, std::string tag, std::string value)
{
    if (tag.compare("RES") == 0)
    {
        setup.resultFolder = value;
        return;
    }
    if (tag.compare("RH") == 0)
    {
        setup.envRelativeHumdity = std::stod(value);
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
    throw Exception(__PRETTY_FUNCTION__, "Unknown tag: \"" + tag + "\". Use -H for help.");
}


void PrintHelp()
{
    std::cout << "Structure for input values: tag=val" << std::endl;
    std::cout << "Possible tags:" << std::endl;
    std::cout << "  res   : result folder" << std::endl;
    std::cout << "  rh    : environmental relative humidity [0.0, 1.0]" << std::endl;
    std::cout << "  ssc   : Stress based shrinkage contribution [0.0, 1.0]" << std::endl;
    std::cout << "  tdry  : drying time [0.0, inf]" << std::endl;
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
    std::cout << "Environmental relative humidity: " << setup.envRelativeHumdity << std::endl;
    std::cout << "Drying time: " << setup.t_dry << std::endl;
    std::cout << "Shrinkage ratio (stress based contribution): " << setup.shrinkageRatio << std::endl;
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

    // Set nodal values ---------------------------------------------------------------------------

    std::cout << "Set nodal values..." << std::endl;

    MPS.SetNodalValues(dofRelativeHumidity, {1.0});
    MPS.SetNodalValues(dofWaterVolumeFraction, {0.1});
    std::cout << "Nodal values set..." << std::endl;

    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsTotal = MPS.CreateCells(groupVolumeElementsTotal, integrationTetrahedron3);
    auto groupSurfaceCellsTotal = MPS.CreateCells(groupSurfaceElementsTotal, integrationTriangle5);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;


    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, MTCoefficientConstant{1.},
                              MTCoefficientConstant{0.01}, MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, 1.,
                              0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});


    auto& integrandMoistureTransportBoundary = MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundary, double>{
            MoistureTransportBoundary{integrandMoistureTransport.GetIntegrand(), dofWaterVolumeFraction,
                                      dofRelativeHumidity, 2., 2., setup.envRelativeHumdity},
            &MoistureTransportBoundary::Gradient, &MoistureTransportBoundary::Stiffness});


    Shrinkage<3> lawShrinkage(dofDisplacements, dofWaterVolumeFraction, dofRelativeHumidity, setup.shrinkageRatio,
                              0.000001);

    auto material = Material::DefaultConcrete();
    material.c = 0.25;
    material.gf *= 0.05;
    GradientDamage<3, Shrinkage<3>> integrandVolume{dofDisplacements, dofNonLocal, material, lawShrinkage};
    int nIp = integrationTetrahedron3.GetNumIntegrationPoints();
    integrandVolume.mKappas = Eigen::MatrixXd::Zero(groupVolumeCellsTotal.Size(), nIp);


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
    double dt_ini = 0.075;
    double dt = dt_ini;


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
    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);


    // Solve --------------------------------------------------------------------------------------


    // Drying Process ----------------------------------------------------

    SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t, setup);

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
            if (dt < dt_ini * 1e-6)
                throw e;

            MPS.MergeDofs(d_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
            MPS.MergeDofs(v_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
            MPS.MergeDofs(a_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);
            continue;
        }

        ++num_converges;
        if (dt < setup.t_dry / 20 && num_converges > 3)
        {
            std::cout << "Increasing timestep" << std::endl;
            dt *= 1.5;
            if (dt > setup.t_dry / 20)
                dt = setup.t_dry / 20;
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
        SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t, setup);
    }


    // Compression Test ----------------------------------------------

    Constraints constraintsCompression;

    auto groupNodesBottom = mesh.NodesAtAxis(eDirection::Z, dofDisplacements);
    auto groupNodesTop = mesh.NodesAtAxis(eDirection::Z, dofDisplacements, 300);
    for (NodeSimple& node : groupNodesBottom)
    {
        Eigen::VectorXd disp = node.GetValues(0);
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{1, 0., 0.}, disp[0]));
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{0., 1, 0.}, disp[1]));
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{0., 0., 1}, disp[2]));
    }
    for (NodeSimple& node : groupNodesTop)
    {
        Eigen::VectorXd disp = node.GetValues(0);
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{1, 0., 0.}, disp[0]));
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{0., 1, 0.}, disp[1]));
        constraintsCompression.Add(dofDisplacements, Direction(node, Eigen::Vector3d{0., 0., 1}, [disp](double t) {
                                       return disp[2] + TimeDependentDisplacement(t);
                                   }));
    }


    t = 0;
    double t_compr = 40;
    dt_ini = t_compr / 400;
    dt = dt_ini;

    MPS.RenumberDofs();

    QSS.SetGlobalTime(t);
    QSS.SetConstraints(constraintsCompression);


    next_plot = dt_ini;
    num_converges = 0;


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
            dt *= 0.5;
            num_converges = 0;
            if (dt < dt_ini * 1e-6)
                throw e;

            MPS.MergeDofs(d_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 0);
            MPS.MergeDofs(v_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 1);
            MPS.MergeDofs(a_MT, {dofWaterVolumeFraction, dofRelativeHumidity}, 2);
            continue;
        }

        ++num_converges;
        if (dt < dt_ini && num_converges > 5)
        {
            std::cout << "Increasing timestep" << std::endl;
            dt *= 1.5;
            if (dt > dt_ini)
                dt = dt_ini;
            num_converges = 0;
        }

        if (t >= next_plot)
        {
            std::cout << std::endl << "Visualize..." << std::endl << std::endl;
            pp.Plot(t + setup.t_dry, false);
            next_plot += dt_ini;
        }
        SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t + setup.t_dry, setup);
    }
}
