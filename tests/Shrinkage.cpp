#include <boost/test/unit_test.hpp>

#include "structures/MultiPhysicsStructure.h"

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"


using namespace NuTo;
using namespace std::placeholders;

BOOST_AUTO_TEST_CASE(Construction)
{
    MultiPhysicsStructure ST;
    ST.CreateUnitMesh(40, 40, InterpolationTriangleQuadratic(), InterpolationTriangleQuadratic(),
                      InterpolationTriangleQuadratic(), InterpolationTrussLobatto(2));


    auto dofs = ST.GetDofs();
    DofType dof_Disp = dofs[0];
    DofType dof_WV = dofs[1];
    DofType dof_RH = dofs[2];
    std::vector<DofType> dofs_ME{dof_Disp};
    std::vector<DofType> dofs_MT{dof_WV, dof_RH};


    Eigen::VectorXd d_Disp = ST.ExtractDofs(dofs_ME, 0);
    Eigen::VectorXd v_Disp = ST.ExtractDofs(dofs_ME, 1);
    Eigen::VectorXd a_Disp = ST.ExtractDofs(dofs_ME, 2);
    Eigen::VectorXd d_MT = ST.ExtractDofs(dofs_MT, 0);
    Eigen::VectorXd v_MT = ST.ExtractDofs(dofs_MT, 1);
    Eigen::VectorXd a_MT = ST.ExtractDofs(dofs_MT, 2);

    const int numDofsWV = ST.GetMesh().NodesTotal(dofs[1]).Size() * dofs[1].GetNum();
    const int numDofsRH = ST.GetMesh().NodesTotal(dofs[2]).Size() * dofs[2].GetNum();

    Constraint::Constraints& constraints = ST.GetConstraints();
    const NodeSimple& nodeBottomLeft = ST.GetMesh().NodeAtCoordinate(Eigen::Vector2d(0., 0.), dof_Disp);
    const NodeSimple& nodeBottomRight = ST.GetMesh().NodeAtCoordinate(Eigen::Vector2d(1., 0.), dof_Disp);
    constraints.Add(dof_Disp,
                    {{nodeBottomLeft, 0, [](double) { return 0.; }},
                     {nodeBottomLeft, 1, [](double) { return 0.; }},
                     {nodeBottomRight, 1, [](double t) {
                          if (t > 0)
                              return t * 0.;
                          return 0.0;
                      }}});

    //    auto test = constraints.BuildUnitConstraintMatrix(dof_Disp, d_Disp.rows());
    //    std::cout << test << std::endl;

    for (int i = 0; i < numDofsWV; ++i)
        d_MT[i] = 0.1;

    for (int i = 0; i < numDofsRH; ++i)
        d_MT[i + numDofsWV] = 1.0;


    ST.MergeDofs(d_MT, dofs_MT, 0);
    ST.MergeDofs(v_MT, dofs_MT, 1);
    ST.MergeDofs(a_MT, dofs_MT, 2);

    double t = 0.;
    double delta_t = 0.005;
    double t_final = 0.1;

    Visualize::PostProcess pp("ShrinkageResults");
    pp.DefineVisualizer("Moisture Transport", ST.GetMoistureTransportCells(), Visualize::AverageHandler());
    pp.Add("Moisture Transport", dof_WV);
    pp.Add("Moisture Transport", dof_RH);
    pp.Add("Moisture Transport", dof_Disp);
    pp.Add("Moisture Transport",
           std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::Stress, ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Stress_Total");
    pp.Add("Moisture Transport", std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::StressMechanics,
                                           ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Stress_Mechanics");
    pp.Add("Moisture Transport", std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::StressShrinkage,
                                           ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Stress_Shrinkage");

    pp.Add("Moisture Transport",
           std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::Strain, ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Strain_Total");
    pp.Add("Moisture Transport", std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::StrainMechanics,
                                           ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Strain_Mechanics");
    pp.Add("Moisture Transport", std::bind(&Integrands::MultiPhysicsMomentumBalance<2>::StrainShrinkage,
                                           ST.GetMultiPhysicsLaw(), _1, std::ref(t)),
           "Strain_Shrinkage");

    pp.Plot(0., false);

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    EigenSparseSolver solver("EigenSparseLU");
    Eigen::VectorXd delta_d = Eigen::VectorXd::Zero(d_MT.rows());

    while (t < t_final)
    {
        t += delta_t;
        int iteration = 0;

        Eigen::VectorXd d_it = d_MT + delta_d;
        Eigen::VectorXd v_it = gamma / (delta_t * beta) * delta_d + (1 - gamma / beta) * v_MT +
                               delta_t * (1 - gamma / (2 * beta)) * a_MT;
        Eigen::VectorXd a_it =
                1 / (delta_t * delta_t * beta) * delta_d - 1 / (delta_t * beta) * v_MT - (1 / (2 * beta) - 1) * a_MT;

        ST.MergeDofs(d_it, dofs_MT, 0);
        ST.MergeDofs(v_it, dofs_MT, 1);
        ST.MergeDofs(a_it, dofs_MT, 2);


        Eigen::VectorXd grad = ToEigen(ST.GradientMoistureTransport(), {dofs[1], dofs[2]}) +
                               ToEigen(ST.GradientMoistureTransportBoundary(), {dofs[1], dofs[2]});

        while (grad.lpNorm<Eigen::Infinity>() > 10e-16)
        {

            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K = ToEigen(ST.StiffnessMoistureTransport(), {dofs[1], dofs[2]}) +
                                            ToEigen(ST.StiffnessMoistureTransportBoundary(), {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> D = ToEigen(ST.DampingMoistureTransport(), {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> H = gamma / (delta_t * beta) * D + K;

            delta_d = solver.Solve(H, -grad);

            d_it += delta_d;
            v_it += gamma / (delta_t * beta) * delta_d;
            a_it += 1 / (delta_t * delta_t * beta) * delta_d;

            ST.MergeDofs(d_it, dofs_MT, 0);
            ST.MergeDofs(v_it, dofs_MT, 1);
            ST.MergeDofs(a_it, dofs_MT, 2);

            grad = ToEigen(ST.GradientMoistureTransport(), {dofs[1], dofs[2]}) +
                   ToEigen(ST.GradientMoistureTransportBoundary(), {dofs[1], dofs[2]});
            std::cout << "Time: " << t << std::endl
                      << "Iteration: " << iteration << std::endl
                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
                      << std::endl;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;


        iteration = 0;
        Eigen::VectorXd B = constraints.GetSparseGlobalRhs(dof_Disp, d_Disp.rows(), t);
        Eigen::VectorXd delta_B = B - constraints.GetSparseGlobalRhs(dof_Disp, d_Disp.rows(), t - delta_t);
        Eigen::SparseMatrix<double> C = constraints.BuildUnitConstraintMatrix(dof_Disp, d_Disp.rows());
        Eigen::SparseMatrix<double> K_Disp = ToEigen(ST.StiffnessMechanics(), {dof_Disp});
        Eigen::VectorXd d_it_disp = d_Disp + delta_B;
        ST.MergeDofs(d_it_disp, {dof_Disp}, 0);

        Eigen::VectorXd grad_Disp =
                ToEigen(ST.GradientMechanics(), {dof_Disp}); // + ToEigen(ST.GradientShrinkage(), {dof_Disp});
        Eigen::VectorXd res_Disp = C.transpose() * grad_Disp;


        while (res_Disp.lpNorm<Eigen::Infinity>() > 10e-4)
        {
            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__,
                                "No convergence - mechanics - residual: " +
                                        std::to_string(res_Disp.lpNorm<Eigen::Infinity>()));

            Eigen::SparseMatrix<double> K_mod = C.transpose() * K_Disp * C;

            Eigen::VectorXd delta_disp = solver.Solve(K_mod, res_Disp);

            d_it_disp += -C * delta_disp;
            ST.MergeDofs(d_it_disp, {dof_Disp}, 0);

            grad_Disp = ToEigen(ST.GradientMechanics(), {dof_Disp}); // + ToEigen(ST.GradientShrinkage(), {dof_Disp});
            res_Disp = C.transpose() * grad_Disp;
        }
        d_Disp = d_it_disp;
        pp.Plot(t, false);
    }
}
