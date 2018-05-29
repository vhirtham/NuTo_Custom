#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
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


using MoistureTransportBoundaryIntegrand =
        Integrands::MoistureTransportBoundary<2, MTCConst<1>, MTCConst<1, 100>, MTCConst<1, 10>>;

int main(int argc, char* argv[])
{
    auto meshGmsh = MeshGmsh("ARAMIS3d.msh");
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};
    std::cout << "Mesh loaded..." << std::endl;

    // Get element groups
    const auto& groupSurfaceElementsBoundary = meshGmsh.GetPhysicalGroup("boundary");
    const auto& groupVolumeElementsMatrix = meshGmsh.GetPhysicalGroup("matrix");
    const auto& groupVolumeElementsGranite = meshGmsh.GetPhysicalGroup("granite");
    const auto& groupVolumeElementsTotal = Unite(groupVolumeElementsMatrix, groupVolumeElementsGranite);

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("Displacements", 3);
    const auto& dofRelativeHumidity = MPS.AddDofType("Relative Humidity", 1);
    const auto& dofWaterVolumeFraction = MPS.AddDofType("Water Volume Fraction", 1);
    std::vector<DofType> dofs_MT{dofWaterVolumeFraction, dofRelativeHumidity};

    // Create interpolations
    auto& interpolationTetrahedronLinear = mesh.CreateInterpolation(InterpolationTetrahedronLinear());
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    auto& interpolationTriangleLinear = mesh.CreateInterpolation(InterpolationTriangleLinear());
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTetrahedron1 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(1));
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));
    const auto& integrationTriangle3 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));

    // Add interpolatins to dof types
    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsMatrix, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsMatrix, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupSurfaceElementsBoundary, interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupSurfaceElementsBoundary, interpolationTriangleQuadratic);
    std::cout << "Dofs created..." << std::endl;

    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);

    // Renumber dofs
    MPS.RenumberDofs();

    // Set nodal values
    MPS.SetNodalValues(dofRelativeHumidity, {1.0});
    MPS.SetNodalValues(dofWaterVolumeFraction, {0.1});
    std::cout << "Nodal values set..." << std::endl;

    // Create Cells
    auto groupVolumeCellsGranite = MPS.CreateCells(groupVolumeElementsGranite, integrationTetrahedron1);
    auto groupVolumeCellsMatrix = MPS.CreateCells(groupVolumeElementsMatrix, integrationTetrahedron3);
    auto groupVolumeCellsTotal = Unite(groupVolumeCellsGranite, groupVolumeCellsMatrix);
    auto groupSurfaceCellsBoundary = MPS.CreateCells(groupSurfaceElementsBoundary, integrationTriangle3);


    // Create integrands

    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, MTCoefficientConstant{1.},
                              MTCoefficientConstant{0.01}, MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, 1.,
                              0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});

    auto& integrandMoistureTransportBoundary =
            MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundaryIntegrand, double>{
                    MoistureTransportBoundaryIntegrand{dofWaterVolumeFraction, dofRelativeHumidity},
                    &MoistureTransportBoundaryIntegrand::Gradient, &MoistureTransportBoundaryIntegrand::Stiffness});


    auto& integrandMechanicsMatrix =
            MPS.AddIntegrand(IntegrandWrapper<MomentumBalanceAdditiveInputOutput<3, Shrinkage<3>>, double>{
                    MomentumBalanceAdditiveInputOutput<3, Shrinkage<3>>{
                            dofDisplacements, Laws::LinearElastic<3>{30000, 0.2},
                            Shrinkage<3>{dofDisplacements, dofWaterVolumeFraction, dofRelativeHumidity, 0.5}},
                    &MomentumBalanceAdditiveInputOutput<3, Shrinkage<3>>::Gradient,
                    &MomentumBalanceAdditiveInputOutput<3, Shrinkage<3>>::Hessian0});

    std::cout << "Integrants created..." << std::endl;

    //    // start test

    //    auto test = MPS.AssembleResidual(dofs_MT,
    //                                     {{&groupVolumeCellsMatrix, &integrandMoistureTransport},
    //                                      {&groupSurfaceCellsBoundary, &integrandMoistureTransportBoundary}});

    //    auto test2 = MPS.AssembleStiffness(dofs_MT,
    //                                       {{&groupVolumeCellsMatrix, &integrandMoistureTransport},
    //                                        {&groupSurfaceCellsBoundary, &integrandMoistureTransportBoundary}});

    //    auto test3 = MPS.AssembleDamping(dofs_MT, {{&groupVolumeCellsMatrix, &integrandMoistureTransport}});

    //    // end test

    Visualize::PostProcess pp("MultiPhysicsResults");
    pp.DefineVisualizer("Matrix", groupVolumeCellsMatrix, Visualize::AverageHandler());
    pp.Add("Matrix", dofDisplacements);
    pp.Add("Matrix", dofRelativeHumidity);
    pp.Add("Matrix", dofWaterVolumeFraction);


    // Build and solve system

    std::cout << "Start time integration scheme..." << std::endl;

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;
    double t = 0.;
    double t_final = 40.;
    double delta_t = t_final / 40.;


    EigenSparseSolver solver("MumpsLU");
    Eigen::VectorXd d_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 0);
    Eigen::VectorXd v_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 1);
    Eigen::VectorXd a_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 2);
    Eigen::VectorXd delta_d_MT = Eigen::VectorXd::Zero(d_MT.rows());

    pp.Plot(t, false);
    while (t < t_final)
    {
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

            std::cout << "Solve system..." << std::endl;
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


        //        iteration = 0;
        //        Eigen::VectorXd B = constraints.GetSparseGlobalRhs(dof_Disp, d_Disp.rows(), t);
        //        Eigen::VectorXd delta_B = B - constraints.GetSparseGlobalRhs(dof_Disp, d_Disp.rows(), t -
        //        delta_t);
        //        Eigen::SparseMatrix<double> C = constraints.BuildUnitConstraintMatrix(dof_Disp, d_Disp.rows());
        //        Eigen::SparseMatrix<double> K_Disp = ToEigen(ST.StiffnessMechanics(), {dof_Disp});
        //        Eigen::VectorXd d_it_disp = d_Disp + delta_B;
        //        ST.MergeDofs(d_it_disp, {dof_Disp}, 0);

        //        Eigen::VectorXd grad_Disp =
        //                ToEigen(ST.GradientMechanics(), {dof_Disp}); // + ToEigen(ST.GradientShrinkage(),
        //                {dof_Disp});
        //        Eigen::VectorXd res_Disp = C.transpose() * grad_Disp;


        //        while (res_Disp.lpNorm<Eigen::Infinity>() > 10e-4)
        //        {
        //            ++iteration;
        //            if (iteration > 20)
        //                throw Exception(__PRETTY_FUNCTION__,
        //                                "No convergence - mechanics - residual: " +
        //                                        std::to_string(res_Disp.lpNorm<Eigen::Infinity>()));

        //            Eigen::SparseMatrix<double> K_mod = C.transpose() * K_Disp * C;

        //            Eigen::VectorXd delta_disp = solver.Solve(K_mod, res_Disp);

        //            d_it_disp += -C * delta_disp;
        //            ST.MergeDofs(d_it_disp, {dof_Disp}, 0);

        //            grad_Disp = ToEigen(ST.GradientMechanics(), {dof_Disp}); // + ToEigen(ST.GradientShrinkage(),
        //            {dof_Disp});
        //            res_Disp = C.transpose() * grad_Disp;
        //        }
        //        d_Disp = d_it_disp;
        std::cout << "Visualize..." << std::endl;
        pp.Plot(t, false);
        std::cout << "Visualization done..." << std::endl;
    }
}
