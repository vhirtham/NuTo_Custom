#include <boost/test/unit_test.hpp>

#include "structures/MultiPhysicsStructure.h"

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/dofs/GlobalDofVector.h"
#include "nuto/mechanics/dofs/GlobalDofMatrixSparse.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"


using namespace NuTo;
// using namespace NuTo::Integrands;


BOOST_AUTO_TEST_CASE(Construction)
{
    MultiPhysicsStructure ST;
    ST.CreateUnitMesh(20, 20, InterpolationTriangleQuadratic(), InterpolationTriangleQuadratic(),
                      InterpolationTriangleQuadratic(), InterpolationTrussLobatto(2));


    auto dofs = ST.GetDofs();
    DofType dof_Disp = dofs[0];
    DofType dof_WV = dofs[1];
    DofType dof_RH = dofs[2];
    std::vector<DofType> dofs_MT{dof_WV, dof_RH};


    Eigen::VectorXd d_MT = ST.ExtractDofs(dofs_MT, 0);
    Eigen::VectorXd v_MT = ST.ExtractDofs(dofs_MT, 1);
    Eigen::VectorXd a_MT = ST.ExtractDofs(dofs_MT, 2);

    const int numDofsWV = ST.GetMesh().NodesTotal(dofs[1]).Size() * dofs[1].GetNum();
    const int numDofsRH = ST.GetMesh().NodesTotal(dofs[2]).Size() * dofs[2].GetNum();


    for (int i = 0; i < numDofsWV; ++i)
    {
        d_MT[i] = 0.1;
    }
    for (int i = 0; i < numDofsRH; ++i)
    {
        d_MT[i + numDofsWV] = 1.0;
    }

    ST.MergeDofs(d_MT, dofs_MT, 0);
    ST.MergeDofs(v_MT, dofs_MT, 1);
    ST.MergeDofs(a_MT, dofs_MT, 2);

    Visualize::PostProcess pp("ShrinkageResults");
    pp.DefineVisualizer("Moisture Transport", ST.GetMoistureTransportCells(), Visualize::AverageHandler());
    pp.Add("Moisture Transport", dof_WV);
    pp.Add("Moisture Transport", dof_RH);
    pp.Plot(0., false);

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    EigenSparseSolver solver("EigenSparseLU");
    Eigen::VectorXd delta_d = Eigen::VectorXd::Zero(d_MT.rows());

    double t = 0.;
    double delta_t = 0.005;
    double t_final = 0.1;
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


        Eigen::VectorXd grad = ToEigen(ST.GradientMoistureTransport().J, {dofs[1], dofs[2]}) +
                               ToEigen(ST.GradientMoistureTransportBoundary().J, {dofs[1], dofs[2]});

        while (grad.lpNorm<Eigen::Infinity>() > 10e-12)
        {

            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K = ToEigen(ST.StiffnessMoistureTransport().JJ, {dofs[1], dofs[2]}) +
                                            ToEigen(ST.StiffnessMoistureTransportBoundary().JJ, {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> D = ToEigen(ST.DampingMoistureTransport().JJ, {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> H = gamma / (delta_t * beta) * D + K;

            delta_d = solver.Solve(H, -grad);

            d_it += delta_d;
            v_it += gamma / (delta_t * beta) * delta_d;
            a_it += 1 / (delta_t * delta_t * beta) * delta_d;

            ST.MergeDofs(d_it, dofs_MT, 0);
            ST.MergeDofs(v_it, dofs_MT, 1);
            ST.MergeDofs(a_it, dofs_MT, 2);

            grad = ToEigen(ST.GradientMoistureTransport().J, {dofs[1], dofs[2]}) +
                   ToEigen(ST.GradientMoistureTransportBoundary().J, {dofs[1], dofs[2]});
            std::cout << "Time: " << t << std::endl
                      << "Iteration: " << iteration << std::endl
                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
                      << std::endl;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;
        pp.Plot(t, false);
    }
}
