#include <boost/test/unit_test.hpp>


#include <iostream>

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"


// NuToCustom includes
#include "MoistureTransportTest.h"
#include "integrands/MoistureTransportCoefficients.h"


//! @brief Checks if the submatrices and subvectors for the two phase system are sized correctly for different
//! interpolation types.
BOOST_AUTO_TEST_CASE(InterpolationCombinations)
{
    for (int i = 1; i < 4; ++i)
        for (int j = 1; j < 4; ++j)
        {
            MoistureTransportTest MTT;
            MTT.CreateUnitMesh(1, InterpolationTrussLobatto(i), InterpolationTrussLobatto(j));
            MTT.Gradient();
            MTT.Stiffness();
            MTT.Damping();
        }
}


////! @brief This test checks if an exception is thrown when the water volume fraction equals the pore volume fraction.
////! This is important because in this case there is no volume left for the gas phase which will lead to undefined
////! behaviour of the gas phase transport equation.
// BOOST_AUTO_TEST_CASE(PoresCompletlyFilledWithWater)
//{

//        MoistureTransportTest<1, MTCConst<1>, MTCConst<0>, MTCConst<1>, MTCConst<2, 10>> MTT(1, 1, 0.2);
//        MTT.CreateUnitMesh(1, InterpolationTrussLobatto(1), InterpolationTrussLobatto(1));
//        BOOST_CHECK_THROW(MTT.Gradient(), Exception);
//        BOOST_CHECK_THROW(MTT.Stiffness(), Exception);
//        BOOST_CHECK_THROW(MTT.Damping(), Exception);
//}
