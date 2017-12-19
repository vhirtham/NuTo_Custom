#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/DirectionEnum.h"
#include "visualize/VisualizeEnum.h"

#include <string>

using namespace NuTo;
using namespace NuTo::Constitutive;
using namespace NuTo::Interpolation;
using namespace NuTo::Node;


// Constants
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// time
constexpr double simulationTime = 20.;
constexpr int numTimesteps = 40;
constexpr int numWriteSteps = 40;

// mechanics
constexpr double youngsModulus = 30000;
constexpr double poissonRatio = 0.2;
constexpr double tensileStrength = 4.;

// mesh
constexpr double cylinderDiameter = 1.;
constexpr double cylinderHeight = 3.;

// iteration
constexpr double tolDisp = 1e-6;
constexpr double tolNLES = 1e-6;
constexpr double tolWF = 1e-10;
constexpr double tolRH = 1e-10;
constexpr int maxIter = 10;


// Simulation
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void DS(std::string simulationName)
{
    Structure S(3);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);
    S.SetVerboseLevel(0);

    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eigen::Vector3d start{0., 0., 0};
    Eigen::Vector3d end{1., 1., 1};
    Eigen::Vector3i subdivisions{8, 8, 30};
    int elementGroup_Id = -1;
    int interpolationType_Id = -1;
    auto meshData = std::tie(elementGroup_Id, interpolationType_Id);
    meshData = MeshGenerator::Grid(S, start, end, subdivisions, eShapeType::BRICK3D);

    assert(elementGroup_Id > -1);
    assert(interpolationType_Id > -1);

    // map unit mesh to cylindrical mesh
    auto cylinderMappingF = MeshGenerator::GetCylinderMapping(cylinderDiameter, cylinderHeight);
    std::vector<int> allNodes_Ids;
    S.NodeGroupGetMembers(S.GroupGetNodesTotal(), allNodes_Ids);
    for (auto& node_Id : allNodes_Ids)
    {
        auto node_Ptr = S.NodeGetNodePtr(node_Id);
        Eigen::VectorXd oldCoordinates = node_Ptr->Get(eDof::COORDINATES);
        Eigen::VectorXd newCoordinates = cylinderMappingF(oldCoordinates);
        node_Ptr->Set(eDof::COORDINATES, newCoordinates);
    }

    // Constitutive laws
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // aggregates - mechanics
    int lawGD_Id = S.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::DENSITY, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::NONLOCAL_RADIUS, 0.2 * 1e-4);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, tensileStrength * 10);
    S.ConstitutiveLawSetDamageLaw(
            lawGD_Id, DamageLawExponential::Create(tensileStrength / youngsModulus, tensileStrength / 0.021));
    S.ElementGroupSetConstitutiveLaw(elementGroup_Id, lawGD_Id);


    // Integration type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationType_Id, eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);

    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);
    //    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::RELATIVEHUMIDITY,
    //                           Interpolation::eTypeOrder::EQUIDISTANT1);
    //    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::WATERVOLUMEFRACTION,
    //                           Interpolation::eTypeOrder::EQUIDISTANT1);


    S.ElementTotalConvertToInterpolationType();


    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    auto& bottomNodesGroup =
            S.GroupGetNodeCoordinateRange(eDirection::Z, -cylinderHeight / 2. - 1e-6, -cylinderHeight / 2. + 1e-6);
    auto& topNodesGroup =
            S.GroupGetNodeCoordinateRange(eDirection::Z, cylinderHeight / 2. - 1e-6, cylinderHeight / 2. + 1e+6);
    assert(bottomNodesGroup.GetNumMembers() > 0);
    assert(topNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Z}));
    S.Constraints().Add(
            Node::eDof::DISPLACEMENTS,
            Constraint::Component(topNodesGroup, {eDirection::Z}, [](double t) { return t / simulationTime * -0.01; }));
    auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(3), 1e-6);
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::X}));
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::Y}));


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int visualizeGroup = elementGroup_Id;
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DAMAGE);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::LOCAL_EQ_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    //    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
    //    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);

    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.SetNumProcessors(4);
    S.CalculateMaximumIndependentSets();

    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, tolDisp);
    NM.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, tolNLES);
    //    NM.SetToleranceResidual(Node::eDof::RELATIVEHUMIDITY, tolRH);
    //    NM.SetToleranceResidual(Node::eDof::WATERVOLUMEFRACTION, tolWF);

    //    NM.AddCalculationStep({Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION});
    //    NM.AddCalculationStep({Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN});
    NM.SetAutomaticTimeStepping(true);
    NM.SetMaxTimeStep(simulationTime / numTimesteps);
    NM.SetMinTimeStep(simulationTime / numTimesteps / 100);
    NM.SetTimeStep(simulationTime / numTimesteps);
    NM.SetPerformLineSearch(false);
    NM.SetMaxNumIterations(maxIter);
    NM.PostProcessing().SetMinTimeStepPlot(simulationTime / numWriteSteps);
    NM.PostProcessing().SetResultDirectory(simulationName, true);
    NM.Solve(simulationTime);
}


// Main
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main()
{
    DS("DamageShrinkage3d");
    return 0;
}
