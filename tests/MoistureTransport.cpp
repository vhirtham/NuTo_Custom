#include <boost/test/unit_test.hpp>

//#define PLOT_RESULTS

#include <iostream>
#include <thread>

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"


// NuToCustom includes
#include "MoistureTransportTest.h"
#include "integrands/MoistureTransportCoefficients.h"

#ifdef PLOT_RESULTS
#include "tools/GNUPlot.h"
#include <chrono>
using namespace std::chrono_literals;
#endif

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void IntegrationTest(int numElements, double delta_t, double t_final, double lWV, double rWV, double lRH, double rRH)
{
    MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq> MTT(1., 0.1, 0.25);

    MTT.CreateUnitMesh(numElements, InterpolationTrussLobatto(2), InterpolationTrussLobatto(2));

    Eigen::VectorXd d, v, a;
    std::tie(d, v, a) = MTT.ExtractDofs();


    auto dofs = MTT.GetDofs();
    const int numDofsWV = MTT.GetMesh().NodesTotal(dofs[0]).Size() * dofs[0].GetNum();
    const int numDofsRH = MTT.GetMesh().NodesTotal(dofs[1]).Size() * dofs[1].GetNum();

    Eigen::VectorXd xWV(numDofsWV), yWV(numDofsWV), xRH(numDofsRH), yRH(numDofsRH);

    for (int i = 0; i < numDofsWV; ++i)
    {
        xWV[i] = static_cast<double>(i) / static_cast<double>(numDofsWV - 1);
        d[i] = lWV + (rWV - lWV) / (numDofsWV - 1) * i;
    }
    for (int i = 0; i < numDofsRH; ++i)
    {
        xRH[i] = static_cast<double>(i) / static_cast<double>(numDofsRH - 1);
        d[i + numDofsWV] = lRH + (rRH - lRH) / (numDofsRH - 1) * i;
    }

    MTT.MergeDofs(d, v, a);

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;

    EigenSparseSolver solver("EigenSparseLU");
    Eigen::VectorXd delta_d = Eigen::VectorXd::Zero(d.rows());
#ifdef PLOT_RESULTS
    GNUPlot plot;
#endif
    double t = 0.;
    while (t < t_final)
    {
        t += delta_t;
        int iteration = 0;

        //        Eigen::VectorXd gradBoundary = ToEigen(MTT.GradientBoundary().J, MTT.GetDofs());
        //        Eigen::MatrixXd stiffBoundary = ToEigen(MTT.StiffnessBoundary().JJ, MTT.GetDofs());
        //        std::cout << stiffBoundary << std::endl;

        Eigen::VectorXd d_it = d + delta_d;
        Eigen::VectorXd v_it =
                gamma / (delta_t * beta) * delta_d + (1 - gamma / beta) * v + delta_t * (1 - gamma / (2 * beta)) * a;
        Eigen::VectorXd a_it =
                1 / (delta_t * delta_t * beta) * delta_d - 1 / (delta_t * beta) * v - (1 / (2 * beta) - 1) * a;
        MTT.MergeDofs(d_it, v_it, a_it);


        Eigen::VectorXd grad =
                ToEigen(MTT.Gradient().J, MTT.GetDofs()) + ToEigen(MTT.GradientBoundary().J, MTT.GetDofs());

        while (grad.lpNorm<Eigen::Infinity>() > 10e-12)
        {

            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K =
                    ToEigen(MTT.Stiffness().JJ, MTT.GetDofs()) + ToEigen(MTT.StiffnessBoundary().JJ, MTT.GetDofs());
            Eigen::SparseMatrix<double> D = ToEigen(MTT.Damping().JJ, MTT.GetDofs());
            Eigen::SparseMatrix<double> H = gamma / (delta_t * beta) * D + K;

            delta_d = solver.Solve(H, -grad);

            d_it += delta_d;
            v_it += gamma / (delta_t * beta) * delta_d;
            a_it += 1 / (delta_t * delta_t * beta) * delta_d;

            MTT.MergeDofs(d_it, v_it, a_it);

            grad = ToEigen(MTT.Gradient().J, MTT.GetDofs()) + ToEigen(MTT.GradientBoundary().J, MTT.GetDofs());
        }
        d = d_it;
        v = v_it;
        a = a_it;
#ifdef PLOT_RESULTS
        // Plotting
        yWV = d.head(numDofsWV);
        yRH = d.tail(numDofsWV);
        plot.Clear();
        plot.AddPlot(std::vector<double>{0, 0}, std::vector<double>{0.0, 1});
        plot.AddPlot(xRH, yRH, {{255, 0, 0}}, eLineType::LINES, "relative humidity");
        plot.AddPlot(xWV, yWV, {{0, 255, 0}}, eLineType::LINES, "water volume fraction");
        plot.Show();
        std::this_thread::sleep_for(20ms);
#endif
    }
}

//! @brief Checks if the submatrices and subvectors for the two phase system are sized correctly for different
//! interpolation types.
BOOST_AUTO_TEST_CASE(InterpolationCombinations)
{
    for (int i = 1; i < 4; ++i)
        for (int j = 1; j < 4; ++j)
        {
            MoistureTransportTest<1, MTCConst<1>, MTCConst<0>, MTCConst<1>, MTCConst<2, 10>> MTT;
            MTT.CreateUnitMesh(1, InterpolationTrussLobatto(i), InterpolationTrussLobatto(j));
            MTT.Gradient();
            MTT.Stiffness();
            MTT.Damping();
        }
}


//! @brief This test checks if an exception is thrown when the water volume fraction equals the pore volume fraction.
//! This is important because in this case there is no volume left for the gas phase which will lead to undefined
//! behaviour of the gas phase transport equation.
BOOST_AUTO_TEST_CASE(PoresCompletlyFilledWithWater)
{
    MoistureTransportTest<1, MTCConst<1>, MTCConst<0>, MTCConst<1>, MTCConst<2, 10>> MTT(1, 1, 0.2);
    MTT.CreateUnitMesh(1, InterpolationTrussLobatto(1), InterpolationTrussLobatto(1));
    BOOST_CHECK_THROW(MTT.Gradient(), Exception);
    BOOST_CHECK_THROW(MTT.Stiffness(), Exception);
    BOOST_CHECK_THROW(MTT.Damping(), Exception);
}

BOOST_AUTO_TEST_CASE(Integrationtest)
{
    //    IntegrationTest<1, MTCConst<4>, MTCConst<2, 10>, MTCConst<20>, MTCConst<1, 10>>(10, 0.0001, 0.05, 0.01, 0.19,
    //    1.0,
    //                                                                                    0.05);
    IntegrationTest<1, MTCConst<4>, MTCConst<2, 100>, MTCConst<20>, MTCConst<1, 10>>(40, 0.0001, 0.05, 0.1, 0.1, 1.0,
                                                                                     1.0);
}
