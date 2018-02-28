#include <boost/test/unit_test.hpp>


// NuToCustom includes
#include "MoistureTransport.h"
#include "integrands/MoistureTransportCoefficients.h"


BOOST_AUTO_TEST_CASE(InterpolationCombinations)
{
    for (int i = 1; i < 4; ++i)
        for (int j = 1; j < 4; ++j)
        {
            MoistureTransportTest<1, MTCConst<1>, MTCConst<1>, MTCConst<1>, MTCConst<2, -1>> MTT;
            MTT.CreateUnitMesh(1, InterpolationTrussLobatto(i), InterpolationTrussLobatto(j));
            MTT.Gradient();
            MTT.Stiffness();
            MTT.Damping();
        }
}
