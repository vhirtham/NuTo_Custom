#include "nuto/math/EigenSparseSolve.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
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

struct TestSetup
{
    double envRelativeHumdity = 0.4;
    double t_dry = 120.; // 400.;
    double initialRH = 1.0;
    double initialWV = 0.18;
    double porosity = 0.2;
    std::string resultFolder{"MTCResults"};
};


void UpdateSetup(TestSetup& setup, std::string tag, std::string value)
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
    throw Exception(__PRETTY_FUNCTION__, "Unknown tag: \"" + tag + "\". Use -H for help.");
}


void PrintHelp()
{
    std::cout << "Structure for input values: tag=val" << std::endl;
    std::cout << "Possible tags:" << std::endl;
    std::cout << "  res   : result folder" << std::endl;
    std::cout << "  rh    : environmental relative humidity [0.0, 1.0]" << std::endl;
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

void PrintSetup(const TestSetup& setup)
{
    std::cout << "-------------" << std::endl;
    std::cout << "Program setup" << std::endl;
    std::cout << "-------------" << std::endl;
    std::cout << "Result folder : " << setup.resultFolder << std::endl;
    std::cout << "Environmental relative humidity: " << setup.envRelativeHumdity << std::endl;
    std::cout << "Drying time: " << setup.t_dry << std::endl;
    std::cout << std::endl << std::endl;
}

TestSetup SetupTestParameters(int argc, char* argv[])
{
    TestSetup setup;

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

    TestSetup setup = SetupTestParameters(argc, argv);
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
    const auto& dofWaterVolumeFraction = MPS.AddDofType("water volume fraction", 1);
    const auto& dofRelativeHumidity = MPS.AddDofType("relative humidity", 1);

    // Create interpolation
    auto& interpolationTetrahedronLinear = mesh.CreateInterpolation(InterpolationTetrahedronLinear());
    auto& interpolationTriangleLinear = mesh.CreateInterpolation(InterpolationTriangleLinear());
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));
    const auto& integrationTriangle5 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupSurfaceElementsTotal, interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupSurfaceElementsTotal, interpolationTriangleQuadratic);

    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);

    // Set nodal values ---------------------------------------------------------------------------

    std::cout << "Set nodal values..." << std::endl;

    MPS.SetNodalValues(dofRelativeHumidity, {setup.initialRH});
    MPS.SetNodalValues(dofWaterVolumeFraction, {setup.initialWV});
    std::cout << "Nodal values set..." << std::endl;

    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsTotal = MPS.CreateCells(groupVolumeElementsTotal, integrationTetrahedron3);
    auto groupSurfaceCellsTotal = MPS.CreateCells(groupSurfaceElementsTotal, integrationTriangle5);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;

    //    use std::abs in coefficient function to avoid negativ coefficients when RH value gets bigger than 1 due to
    //            oscillations

    //    auto diffusionWV = MTCoefficientLinear{0., -10., 0., 0., 11.};
    auto diffusionWV = MTCoefficientQuadratic{{{0., 0.}}, {{-10., 0.}}, {{0., 0.}}, {{0., 0.}}, 10.1};
    auto diffusionRH = MTCoefficientConstant{0.};
    auto massExchange = MTCoefficientConstant{10.};
    auto sorptionIso = MTCoefficientLinear{0.0, setup.initialWV, 0., 0., 0.};

    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, diffusionWV, diffusionRH, massExchange,
                              sorptionIso, 1., 0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});


    auto& integrandMoistureTransportBoundary = MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundary, double>{
            MoistureTransportBoundary{integrandMoistureTransport.GetIntegrand(), dofWaterVolumeFraction,
                                      dofRelativeHumidity, 8., 0., setup.envRelativeHumdity},
            &MoistureTransportBoundary::Gradient, &MoistureTransportBoundary::Stiffness});


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


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    PostProcess pp(setup.resultFolder);
    pp.DefineVisualizer("Volume", groupVolumeCellsTotal, AverageHandler());
    pp.Add("Volume", dofWaterVolumeFraction);
    pp.Add("Volume", dofRelativeHumidity);

    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);


    // Solve --------------------------------------------------------------------------------------


    double delta_plot = dt;
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


        if (std::abs(1. - t) <= 1.e-6)
            dt = 1.;
        if (std::abs(2. - t) <= 1.e-6)
            dt = 2.;
        //        if (dt < setup.t_dry / 20)
        //        {
        //            std::cout << "Increasing timestep" << std::endl;
        //            dt *= 1.5;
        //            if (dt > setup.t_dry / 20)
        //                dt = setup.t_dry / 20;
        //            if (t + dt > setup.t_dry)
        //                dt = setup.t_dry - t;
        //            num_converges = 0;
        //        }

        if (t >= next_plot)
        {
            std::cout << std::endl << "Visualize..." << std::endl << std::endl;
            pp.Plot(t, false);
            next_plot += delta_plot;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;
    }
}
