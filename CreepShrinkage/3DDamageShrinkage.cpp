#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/elements/ElementBase.h"
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
constexpr int numTimesteps = 160;
constexpr int numWriteSteps = 160;

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

void DS(double RH_Air, bool useStressbased, std::string simulationName)
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
    Eigen::Vector3i subdivisions{6, 6, 18};
    int elementGroup_Id = -1;
    int interpolationType_Id = -1;
    auto meshData = std::tie(elementGroup_Id, interpolationType_Id);
    meshData = MeshGenerator::Grid(S, start, end, subdivisions, eShapeType::BRICK3D);

    assert(elementGroup_Id > -1);
    assert(interpolationType_Id > -1);

    std::vector<int> allNodes_Ids;
    S.NodeGroupGetMembers(S.GroupGetNodesTotal(), allNodes_Ids);


    // Constitutive laws
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // mechanics
    //    int lawGD_Id = S.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int lawGD_Id = S.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);

    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);

    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::DENSITY, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::NONLOCAL_RADIUS, 0.2 * 1e-5);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, tensileStrength * 10);
    S.ConstitutiveLawSetDamageLaw(lawGD_Id,
                                  DamageLawExponential::Create(tensileStrength / youngsModulus, tensileStrength / 0.1));

    // moisture transport

    // matrix - moisture transport
    int lawMT_Id = S.ConstitutiveLawCreate(eConstitutiveType::MOISTURE_TRANSPORT);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 1.0e-3);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 1.0e-1);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH, 3.9e-6);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV, 1.17e-1);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_RH, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_WV, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION,
                                        0.56);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION,
                                        0.26);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::MASS_EXCHANGE_RATE, 0.e-8);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::PORE_VOLUME_FRACTION, 0.25);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR, 0.0173);
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMT_Id, eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,
            (Eigen::Vector4d() << 0., 0.19692057340725558, -0.28253538941816925, 0.22661481601091368).finished());
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMT_Id, eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,
            (Eigen::Vector4d() << 0., 0.26719233184420238, -0.41030868184510738, 0.32511635000090505).finished());

    // Strain based shrinkage
    int lawSSaID = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    S.ConstitutiveLawSetParameterDouble(lawSSaID, eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,
                                        youngsModulus * 4. / 3. * 10e6);
    S.ConstitutiveLawSetParameterDouble(lawSSaID, eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                        youngsModulus * 2. / 3. * 10e6);
    S.ConstitutiveLawSetParameterDouble(lawSSaID, eConstitutiveParameter::TEMPERATURE, 293);

    // Stress based shrinkage
    int lawSSeID = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED);
    S.ConstitutiveLawSetParameterDouble(lawSSeID, eConstitutiveParameter::TEMPERATURE, 293);

    // Additive laws
    int lawAO_Id = S.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);
    auto lawAO_Ptr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAO_Id));
    int lawAIE_Id = S.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    auto lawAIE_Ptr = static_cast<NuTo::AdditiveInputExplicit*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_Id));
    if (useStressbased)
    {
        lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawGD_Id));
        lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawSSeID));
        lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMT_Id));
    }
    else
    {
        lawAIE_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawGD_Id));
        lawAIE_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawSSaID), eInput::ENGINEERING_STRAIN);
        lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_Id));
        lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMT_Id));
    }

    S.ElementGroupSetConstitutiveLaw(elementGroup_Id, lawAO_Id);


    // Integration type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationType_Id, eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);

    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::RELATIVEHUMIDITY, Interpolation::eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationType_Id, Node::eDof::WATERVOLUMEFRACTION,
                           Interpolation::eTypeOrder::EQUIDISTANT2);


    S.ElementTotalConvertToInterpolationType();

    // Create necessary groups before mesh mapping
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto& bottomNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::Z, 0. - 1e-6, 0. + 1e-6);
    auto& topNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::Z, 1. - 1e-6, 1. + 1e+6);


    // map unit mesh to cylindrical mesh and warp to avoid biforcation
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto cylinderMappingF = MeshGenerator::GetCylinderMapping(cylinderDiameter, cylinderHeight);
    auto antiBiforkationMapping = [](Eigen::VectorXd oldCoordinates) -> Eigen::VectorXd {
        oldCoordinates[2] = oldCoordinates[2] +
                            (((oldCoordinates[0] / cylinderDiameter) * (oldCoordinates[1] / cylinderDiameter) +
                              (oldCoordinates[1] / cylinderDiameter)) *
                             ((oldCoordinates[2] / cylinderHeight) + 0.5)) *
                                    0.002;
        return oldCoordinates;
    };

    for (auto& node_Id : allNodes_Ids)
    {
        auto node_Ptr = S.NodeGetNodePtr(node_Id);
        Eigen::VectorXd oldCoordinates = node_Ptr->Get(Node::eDof::COORDINATES);
        Eigen::VectorXd newCoordinates = antiBiforkationMapping(cylinderMappingF(oldCoordinates));
        node_Ptr->Set(Node::eDof::COORDINATES, newCoordinates);
    }

    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double initialRH = 0.95;
    int BoundaryNodesGroup = S.GroupCreateNodeGroup();
    S.GroupAddNodeFunction(BoundaryNodesGroup, [](NodeBase* node) -> bool {
        auto coords = node->Get(Node::eDof::COORDINATES);
        auto r = std::sqrt(coords[0] * coords[0] + coords[1] * coords[1]);
        if (r > std::abs(cylinderDiameter / 2.) - 1e-6)
            return true;
        return false;
    });
    int ElementsAtBoundariesGroup = S.GroupCreateElementGroup();
    S.GroupAddElementsFromNodes(ElementsAtBoundariesGroup, BoundaryNodesGroup, false);


    S.GroupAddElementsFromNodes(ElementsAtBoundariesGroup, BoundaryNodesGroup, false);

    auto boundaryControlNodePtr = S.NodeGetNodePtr(S.NodeCreateDOFs({Node::eDof::RELATIVEHUMIDITY}));

    int boundaryElementsGroup =
            S.BoundaryElementsCreate(ElementsAtBoundariesGroup, BoundaryNodesGroup, boundaryControlNodePtr);

    S.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                        NuTo::Constraint::Value(*boundaryControlNodePtr, [initialRH, RH_Air](double rTime) -> double {
                            if (rTime == 0.) // || rTime > SIMULATIONTIME / 2.0)
                                return initialRH;
                            else if (rTime < simulationTime / 20.)
                                return initialRH - rTime / (simulationTime / 20.) * (initialRH - RH_Air);
                            else
                                return RH_Air;
                        }));

    std::vector<int> boundaryElementIDsVector;
    S.ElementGroupGetMembers(boundaryElementsGroup, boundaryElementIDsVector);


    for (int elementId : boundaryElementIDsVector)
    {
        S.ElementGetElementPtr(elementId)->SetIntegrationType(
                *S.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType2D4NGauss9Ip));
    }

    // Initial Values
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    std::vector<int> nodeIDVectorMatrix;

    int allElementsGroup_Id = S.GroupCreateElementGroup();
    S.GroupAddElementsTotal(allElementsGroup_Id);

    S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(allElementsGroup_Id), nodeIDVectorMatrix);
    double initialWV = S.ConstitutiveLawGetEquilibriumWaterVolumeFraction(
            lawMT_Id, initialRH,
            S.ConstitutiveLawGetParameterFullVectorDouble(
                    lawMT_Id, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
    for (auto nodeID : nodeIDVectorMatrix)
    {
        auto nodePtr = S.NodeGetNodePtr(nodeID);
        if (nodePtr->IsDof(NuTo::Node::eDof::RELATIVEHUMIDITY))
            nodePtr->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, initialRH);
        if (nodePtr->IsDof(NuTo::Node::eDof::WATERVOLUMEFRACTION))
            nodePtr->Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, initialWV);
    }


    // Initialize Static Data
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto lawMT_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(lawMT_Id);
    std::vector<int> elementIDsMatrixVector;

    S.ElementGroupGetMembers(allElementsGroup_Id, elementIDsMatrixVector);

    for (unsigned int i = 0; i < elementIDsMatrixVector.size(); i++)
    {
        for (int theIP = 0; theIP < S.ElementGetElementPtr(elementIDsMatrixVector[i])->GetNumIntegrationPoints();
             theIP++)
        {
            IPAdditiveOutput& ipLawAO =
                    dynamic_cast<IPAdditiveOutput&>(S.ElementGetElementPtr(i)->GetIPData().GetIPConstitutiveLaw(theIP));

            StaticData::DataMoistureTransport& moistureData =
                    ipLawAO.GetSublawData<MoistureTransport>(lawMT_Ptr).GetData(); // finally the data.
            moistureData.SetLastSorptionCoeff(lawMT_Ptr->GetParameterFullVectorDouble(
                    eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
            moistureData.SetCurrentSorptionCoeff(lawMT_Ptr->GetParameterFullVectorDouble(
                    eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
            moistureData.SetLastRelHumValue(initialRH);
            moistureData.SetDesorption(true);
        }
    }
    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    assert(bottomNodesGroup.GetNumMembers() > 0);
    assert(topNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Z}));
    S.Constraints().Add(
            Node::eDof::DISPLACEMENTS,
            Constraint::Component(topNodesGroup, {eDirection::Z}, [](double t) { return t / simulationTime * -0.02; }));
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
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);

    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.SetNumProcessors(4);
    S.CalculateMaximumIndependentSets();

    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, tolDisp);
    NM.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, tolNLES);
    NM.SetToleranceResidual(Node::eDof::RELATIVEHUMIDITY, tolRH);
    NM.SetToleranceResidual(Node::eDof::WATERVOLUMEFRACTION, tolWF);

    NM.AddCalculationStep({Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION});
    NM.AddCalculationStep({Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN});
    NM.SetAutomaticTimeStepping(true);
    NM.SetMaxTimeStep(simulationTime / numTimesteps);
    NM.SetMinTimeStep(simulationTime / numTimesteps / 1000);
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
    DS(0.3, true, "DamageShrinkage3d");
    return 0;
}
