#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"

#include "nuto/mechanics/cell/SimpleAssembler.h"

#include "structures/MultiPhysicsStructureNew.h"

#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportCoefficients.h"
#include "tools/IntegrandWrapper.h"

using namespace NuTo;
using namespace std::placeholders;

int main(int argc, char* argv[])
{
    auto meshGmsh = MeshGmsh("ARAMIS3d.msh");
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};

    // Get element groups
    const auto& groupVolumeElementsMatrix = meshGmsh.GetPhysicalGroup("matrix");
    const auto& groupVolumeElementsGranite = meshGmsh.GetPhysicalGroup("granite");
    const auto& groupVolumeElementsTotal = mesh.ElementsTotal();

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("Displacements", 3);
    const auto& dofRelativeHumidity = MPS.AddDofType("Relative Humidity", 1);
    const auto& dofWaterVolumeFraction = MPS.AddDofType("Water Volume Fraction", 1);

    // Create interpolations
    auto& interpolationTetrahedronLinear = mesh.CreateInterpolation(InterpolationTetrahedronLinear());
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());

    // Create integrations
    const auto& integrationTetrahedron1 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(1));
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(3));

    // Add interpolatins to dof types
    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsMatrix, interpolationTetrahedronLinear);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsMatrix, interpolationTetrahedronLinear);

    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);

    // Renumber dofs
    MPS.RenumberDofs();

    // Set nodal values
    MPS.SetNodalValues(dofRelativeHumidity, {1.0});
    MPS.SetNodalValues(dofWaterVolumeFraction, {0.1});

    // Create Cells
    auto groupVolumeCellsGranite = MPS.CreateCells(groupVolumeElementsGranite, integrationTetrahedron1);
    auto groupVolumeCellsMatrix = MPS.CreateCells(groupVolumeElementsMatrix, integrationTetrahedron3);
    auto groupVolumeCellsTotal = MPS.GetAllCells();


    // test


    auto& integrand = MPS.AddIntegrand(IntegrandWrapper<Integrands::MoistureTransport>(
            &Integrands::MoistureTransport::Gradient, dofWaterVolumeFraction, dofRelativeHumidity,
            MTCoefficientConstant{1.}, MTCoefficientConstant{0.01}, MTCoefficientConstant{0.},
            MTCoefficientConstant{0.1}, 1., 0.1, 0.2));

    auto test = SimpleAssembler(MPS.GetDofInfo())
                        .BuildVector(groupVolumeCellsMatrix, {dofRelativeHumidity, dofWaterVolumeFraction},
                                     std::bind(&IntegrandWrapperInterface::Residual, &integrand, _1, 0.));

    // end test


    Visualize::PostProcess pp("MultiPhysicsResults");
    pp.DefineVisualizer("Matrix", groupVolumeCellsMatrix, Visualize::AverageHandler());
    pp.Add("Matrix", dofDisplacements);
    pp.Add("Matrix", dofRelativeHumidity);
    pp.Add("Matrix", dofWaterVolumeFraction);
    pp.Plot(0., false);

    int a = 0;
}
