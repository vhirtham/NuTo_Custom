#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"

#include "structures/MultiPhysicsStructureNew.h"
#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportBoundary.h"
#include "integrands/MoistureTransportCoefficients.h"
#include "integrands/MomentumBalanceAdditiveInputOutput.h"
#include "integrands/Shrinkage.h"
#include "tools/IntegrandWrapper.h"

using namespace NuTo;
using namespace NuTo::Integrands;
using namespace std::placeholders;

using IntegrandMechanicsMatrix = MomentumBalanceAdditiveInputOutput<3, Shrinkage<3>>;

void ApplyDisplacementConstraints(MultiPhysicsStructure& MPS, const DofType& dofDisplacements)
{
    auto& constraints = MPS.GetConstraints();
    auto& mesh = MPS.GetMesh();

    const NodeSimple& nodeFrontLeftBottom = mesh.NodeAtCoordinate(Eigen::Vector3d{0., 0., 0.}, dofDisplacements);
    const NodeSimple& nodeFrontLeftTop = mesh.NodeAtCoordinate(Eigen::Vector3d{0., 200., 0.}, dofDisplacements);
    const NodeSimple& nodeFrontRightBottom = mesh.NodeAtCoordinate(Eigen::Vector3d{200., 0., 0.}, dofDisplacements);

    constraints.Add(dofDisplacements,
                    {{nodeFrontLeftBottom, 0, [](double) { return 0.; }},
                     {nodeFrontLeftBottom, 1, [](double) { return 0.; }},
                     {nodeFrontLeftBottom, 2, [](double) { return 0.; }},
                     {nodeFrontLeftTop, 0, [](double) { return 0.; }},
                     {nodeFrontLeftTop, 2, [](double) { return 0.; }},
                     {nodeFrontRightBottom, 1, [](double) { return 0.; }},
                     {nodeFrontRightBottom, 2, [](double) { return 0.; }}});
    MPS.RenumberDofs();
}


int main(int argc, char* argv[])
{
#ifdef NDEBUG
    auto meshGmsh = MeshGmsh("ARAMIS3d.msh");
#else
    auto meshGmsh = MeshGmsh("ARAMIS3d_coarse.msh");
#endif
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};
    std::cout << "Mesh loaded..." << std::endl;

    // Get element groups
    const auto& groupSurfaceElementsLeftMatrixBoundary = meshGmsh.GetPhysicalGroup("leftmatrixboundary");
    const auto& groupSurfaceElementsRightMatrixBoundary = meshGmsh.GetPhysicalGroup("rightmatrixboundary");
    const auto& groupSurfaceElementsTopMatrixBoundary = meshGmsh.GetPhysicalGroup("topmatrixboundary");
    const auto& groupSurfaceElementsBottomMatrixBoundary = meshGmsh.GetPhysicalGroup("bottommatrixboundary");
    const auto& groupSurfaceElementsFrontMatrixBoundary = meshGmsh.GetPhysicalGroup("frontmatrixboundary");
    const auto& groupSurfaceElementsBackMatrixBoundary = meshGmsh.GetPhysicalGroup("backmatrixboundary");
    const auto& groupSurfaceElementsTotalMatrixBoundary =
            Unite(groupSurfaceElementsLeftMatrixBoundary, groupSurfaceElementsRightMatrixBoundary,
                  groupSurfaceElementsTopMatrixBoundary, groupSurfaceElementsBottomMatrixBoundary,
                  groupSurfaceElementsFrontMatrixBoundary, groupSurfaceElementsBackMatrixBoundary);
    const auto& groupVolumeElementsMatrix = meshGmsh.GetPhysicalGroup("matrix");
    const auto& groupVolumeElementsGranite = meshGmsh.GetPhysicalGroup("granite");
    const auto& groupVolumeElementsTotal = Unite(groupVolumeElementsMatrix, groupVolumeElementsGranite);

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("Displacements", 3);
    const auto& dofRelativeHumidity = MPS.AddDofType("Relative Humidity", 1);
    const auto& dofWaterVolumeFraction = MPS.AddDofType("Water Volume Fraction", 1);
    std::vector<DofType> dofs_MT{dofWaterVolumeFraction, dofRelativeHumidity};

    // Create interpolations
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTetrahedron5 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));
    const auto& integrationTriangle5 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsMatrix, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsMatrix, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupSurfaceElementsTotalMatrixBoundary,
                        interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupSurfaceElementsTotalMatrixBoundary,
                        interpolationTriangleQuadratic);

    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);


    // Set nodal values ---------------------------------------------------------------------------

    std::cout << "Set nodal values..." << std::endl;

    MPS.SetNodalValues(dofRelativeHumidity, {1.0});
    MPS.SetNodalValues(dofWaterVolumeFraction, {0.1});
    std::cout << "Nodal values set..." << std::endl;

    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsGranite = MPS.CreateCells(groupVolumeElementsGranite, integrationTetrahedron5);
    auto groupVolumeCellsMatrix = MPS.CreateCells(groupVolumeElementsMatrix, integrationTetrahedron5);
    auto groupVolumeCellsTotal = Unite(groupVolumeCellsGranite, groupVolumeCellsMatrix);
    auto groupSurfaceCellsBoundary = MPS.CreateCells(groupSurfaceElementsTotalMatrixBoundary, integrationTriangle5);


    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;

    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, MTCoefficientConstant{1.},
                              MTCoefficientConstant{0.01}, MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, 1.,
                              0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});


    auto& integrandMoistureTransportBoundary = MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundary, double>{
            MoistureTransportBoundary{integrandMoistureTransport.GetIntegrand(), dofWaterVolumeFraction,
                                      dofRelativeHumidity, 2., 2., 0.4},
            &MoistureTransportBoundary::Gradient, &MoistureTransportBoundary::Stiffness});


    auto& integrandMechanicsMatrix = MPS.AddIntegrand(IntegrandWrapper<IntegrandMechanicsMatrix, double>{
            IntegrandMechanicsMatrix{dofDisplacements, Laws::LinearElastic<3>{30.e9, 0.2},
                                     Shrinkage<3>{dofDisplacements, dofWaterVolumeFraction, dofRelativeHumidity, 0.5}},
            &IntegrandMechanicsMatrix::Gradient, &IntegrandMechanicsMatrix::Hessian0});


    auto lawGranite = Laws::LinearElastic<3>{60.e9, 0.2};
    auto& integrandMechanicsGranite = MPS.AddIntegrand(
            IntegrandWrapper<MomentumBalance<3>, double>{MomentumBalance<3>{dofDisplacements, lawGranite},
                                                         &MomentumBalance<3>::Gradient, &MomentumBalance<3>::Hessian0});


    // Setup time integration ---------------------------------------------------------------------

    std::cout << "Setup time integration scheme..." << std::endl;

    double t = 0.;
    double t_final = 40.;
    double delta_t = t_final / 40.;


    // Newmark

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    // Quasi static

    ApplyDisplacementConstraints(MPS, dofDisplacements);

    TimeDependentProblem TDP{&mesh};
    TDP.AddGradientFunction(groupVolumeCellsMatrix, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMechanicsMatrix.Residual(cipd, dt);
    });
    TDP.AddHessian0Function(groupVolumeCellsMatrix, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMechanicsMatrix.Stiffness(cipd, dt);
    });
    TDP.AddGradientFunction(groupVolumeCellsGranite, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMechanicsGranite.Residual(cipd, dt);
    });
    TDP.AddHessian0Function(groupVolumeCellsGranite, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMechanicsGranite.Stiffness(cipd, dt);
    });
    QuasistaticSolver QSS{TDP, dofDisplacements};
    QSS.SetGlobalTime(t);
    QSS.SetConstraints(MPS.GetConstraints());


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    Visualize::PostProcess pp("MultiPhysicsResults");
    pp.DefineVisualizer("Matrix", groupVolumeCellsMatrix, Visualize::AverageHandler());
    pp.Add("Matrix", dofDisplacements);
    pp.Add("Matrix", dofRelativeHumidity);
    pp.Add("Matrix", dofWaterVolumeFraction);
    pp.Add("Matrix",
           std::bind(&IntegrandMechanicsMatrix::Stress, &(integrandMechanicsMatrix.GetIntegrand()), _1, std::ref(t)),
           "Stress_Total");
    pp.Add("Matrix",
           std::bind(&IntegrandMechanicsMatrix::Strain, &(integrandMechanicsMatrix.GetIntegrand()), _1, std::ref(t)),
           "Strain_Total");
    pp.Add("Matrix", std::bind(&IntegrandMechanicsMatrix::StressMechanics, &(integrandMechanicsMatrix.GetIntegrand()),
                               _1, std::ref(t)),
           "Stress_Mechanics");
    pp.Add("Matrix", std::bind(&IntegrandMechanicsMatrix::StrainMechanics, &(integrandMechanicsMatrix.GetIntegrand()),
                               _1, std::ref(t)),
           "Strain_Mechanics");
    pp.Add("Matrix", std::bind(&IntegrandMechanicsMatrix::StressAdditive, &(integrandMechanicsMatrix.GetIntegrand()),
                               _1, std::ref(t)),
           "Stress_Shrinkage");
    pp.Add("Matrix", std::bind(&IntegrandMechanicsMatrix::StrainAdditive, &(integrandMechanicsMatrix.GetIntegrand()),
                               _1, std::ref(t)),
           "Strain_Shrinkage");


    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);


    // Solve --------------------------------------------------------------------------------------

    std::cout << "Extract nodal values..." << std::endl;

    Eigen::VectorXd d_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 0);
    Eigen::VectorXd v_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 1);
    Eigen::VectorXd a_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 2);
    Eigen::VectorXd delta_d_MT = Eigen::VectorXd::Zero(d_MT.rows());

    std::cout << "Start time integration scheme..." << std::endl;
    EigenSparseSolver solver("MumpsLU");
    while (t < t_final)
    {

        // Moisture Transport - Newmark/Crank-Nicholsen

        t += delta_t;
        int iteration = 0;

        Eigen::VectorXd d_it = d_MT + delta_d_MT;
        Eigen::VectorXd v_it = gamma / (delta_t * beta) * delta_d_MT + (1 - gamma / beta) * v_MT +
                               delta_t * (1 - gamma / (2 * beta)) * a_MT;
        Eigen::VectorXd a_it =
                1 / (delta_t * delta_t * beta) * delta_d_MT - 1 / (delta_t * beta) * v_MT - (1 / (2 * beta) - 1) * a_MT;

        MPS.MergeDofs(d_it, dofs_MT, 0);
        MPS.MergeDofs(v_it, dofs_MT, 1);
        MPS.MergeDofs(a_it, dofs_MT, 2);


        Eigen::VectorXd grad =
                MPS.AssembleResidual(dofs_MT,
                                     {{&groupVolumeCellsMatrix, &integrandMoistureTransport},
                                      {&groupSurfaceCellsBoundary, &integrandMoistureTransportBoundary}});

        while (grad.lpNorm<Eigen::Infinity>() > 10e-13)
        {

            std::cout << "Assemble moisture problem..." << std::endl;
            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K =
                    MPS.AssembleStiffness(dofs_MT,
                                          {{&groupVolumeCellsMatrix, &integrandMoistureTransport},
                                           {&groupSurfaceCellsBoundary, &integrandMoistureTransportBoundary}});

            Eigen::SparseMatrix<double> D =
                    MPS.AssembleDamping(dofs_MT, {{&groupVolumeCellsMatrix, &integrandMoistureTransport}});
            Eigen::SparseMatrix<double> H = gamma / (delta_t * beta) * D + K;

            std::cout << "Solve moisture problem..." << std::endl;
            delta_d_MT = solver.Solve(H, -grad);

            d_it += delta_d_MT;
            v_it += gamma / (delta_t * beta) * delta_d_MT;
            a_it += 1 / (delta_t * delta_t * beta) * delta_d_MT;

            MPS.MergeDofs(d_it, dofs_MT, 0);
            MPS.MergeDofs(v_it, dofs_MT, 1);
            MPS.MergeDofs(a_it, dofs_MT, 2);

            grad = MPS.AssembleResidual(dofs_MT,
                                        {{&groupVolumeCellsMatrix, &integrandMoistureTransport},
                                         {&groupSurfaceCellsBoundary, &integrandMoistureTransportBoundary}});

            std::cout << "Time: " << t << std::endl
                      << "Iteration: " << iteration << std::endl
                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
                      << std::endl;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;


        // Mechanics - Quasi static

        std::cout << "Assemble and solve mechanics problem..." << std::endl;
        QSS.mTolerance = 10e-4;
        QSS.DoStep(t, "MumpsLU");

        std::cout << "Visualize..." << std::endl;
        pp.Plot(t, false);
    }
}
