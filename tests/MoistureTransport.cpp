#include <boost/test/unit_test.hpp>

#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(IntegrationTest)
{
    // Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr int numElements = 100;
    MeshFem mesh = UnitMeshFem::CreateLines(numElements);

    DofType dofRH("relativeHumidity", 1);
    DofType dofWV("waterVolumeFraction", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLobatto(2));
    AddDofInterpolation(&mesh, dofRH, interpolation);
    AddDofInterpolation(&mesh, dofWV, interpolation);

    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    Constraint::Constraints constraints;


    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofNumbering::DofInfo dofRHInfo = DofNumbering::Build(mesh.NodesTotal(dofRH), dofRH, constraints);
    DofNumbering::DofInfo dofWVInfo = DofNumbering::Build(mesh.NodesTotal(dofWV), dofWV, constraints);
}
