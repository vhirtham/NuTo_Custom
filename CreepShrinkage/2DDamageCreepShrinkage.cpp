#include "mechanics/DirectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/ShrinkageCapillaryStrainBased.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "visualize/VisualizeEnum.h"

using namespace Eigen;
using namespace NuTo;
using namespace NuTo::Constitutive;
using namespace NuTo::Constraint;

// Constants
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// time
constexpr double simulationTime = 240.;
constexpr int numTimesteps = 480;
constexpr int numWriteSteps = 60;

// mechanics
constexpr double youngsModulus = 30000;
constexpr double poissonRatio = 0.2;
constexpr double tensileStrength = 4.;


// iteration
constexpr double tolDisp = 1e-6;
constexpr double tolNLES = 1e-6;
constexpr double tolWF = 1e-10;
constexpr double tolRH = 1e-10;
constexpr int maxIter = 10;


// Simulation
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//! @brief Simulation with damage, creep and shrinkage
void DCS(std::string simulationName)
{
    Structure S(2);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);
    S.SetVerboseLevel(0);

    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto meshData = S.ImportFromGmsh("Meso.msh");

    int elementGroupMatrixID = meshData[0].first;
    int interpolationTypeMatrixID = meshData[0].second;

    int elementGroupAggregatesID = meshData[1].first;
    int interpolationTypeAggregatesID = meshData[1].second;

    int allElementsGroup = S.GroupCreateElementGroup();
    S.GroupAddElementsTotal(allElementsGroup);

    // Section
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.ElementTotalSetSection(SectionPlane::Create(1.0, true));

    // Constitutive laws
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix - mechanics
    int lawGD_Id = S.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::DENSITY, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::NONLOCAL_RADIUS, 0.2 * 1e-4);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, tensileStrength * 10);
    S.ConstitutiveLawSetDamageLaw(
            lawGD_Id, DamageLawExponential::Create(tensileStrength / youngsModulus, tensileStrength / 0.021));

    //    int lawGD_Id = S.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    //    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    //    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);


    // matrix - moisture transport
    int lawMT_Id = S.ConstitutiveLawCreate(eConstitutiveType::MOISTURE_TRANSPORT);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 1.0e-3);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 1.0e1);
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
    int lawShr_Id = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    S.ConstitutiveLawSetParameterDouble(lawShr_Id, eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,
                                        youngsModulus * 4. / 3. * 10.e6 / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShr_Id, eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                        youngsModulus * 2. / 3. * 10.e6 / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShr_Id, eConstitutiveParameter::TEMPERATURE, 293);


    // matrix
    int lawAO_Id = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int lawAIE_Id = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    auto lawAO_Ptr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAO_Id));
    auto lawAIE_Ptr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_Id));

    lawAIE_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawGD_Id));
    lawAIE_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawShr_Id), eInput::ENGINEERING_STRAIN);
    lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_Id));
    lawAO_Ptr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMT_Id));


    S.ElementGroupSetConstitutiveLaw(elementGroupMatrixID, lawAO_Id);


    // aggregates - mechanics
    int lawLE_Id = S.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(lawLE_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus / 200. * 400);
    S.ConstitutiveLawSetParameterDouble(lawLE_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    S.ElementGroupSetConstitutiveLaw(elementGroupAggregatesID, lawLE_Id);

    // Integration type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationTypeMatrixID, eIntegrationType::IntegrationType2D4NGauss9Ip);

    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::DISPLACEMENTS,
                           Interpolation::eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::NONLOCALEQSTRAIN,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::RELATIVEHUMIDITY,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::WATERVOLUMEFRACTION,
                           Interpolation::eTypeOrder::EQUIDISTANT1);

    // aggregates
    S.InterpolationTypeAdd(interpolationTypeAggregatesID, Node::eDof::DISPLACEMENTS,
                           Interpolation::eTypeOrder::EQUIDISTANT2);

    S.ElementTotalConvertToInterpolationType();


    // Initial Values
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    std::vector<int> nodeIDVectorMatrix;

    S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(elementGroupMatrixID), nodeIDVectorMatrix);
    double initialRH = 0.95;
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

    S.ElementGroupGetMembers(elementGroupMatrixID, elementIDsMatrixVector);

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
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.2, 0.2);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS,
    //                        Constraint::Component(rightNodesGroup, {eDirection::X},
    //                                              [](double t) { return t / simulationTime * -0.002 * -0.00; }));
    auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(2));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));


    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int upperSideNodesGroup = S.GroupCreateNodeGroup();
    int upperSideElementsGroup = S.GroupCreateElementGroup();

    S.GroupAddNodeCoordinateRange(upperSideNodesGroup, eDirection::Y, 0. - 1.e-6, 0. + 1.e-6);
    S.GroupAddElementsFromNodes(upperSideElementsGroup, upperSideNodesGroup, false);

    auto boundaryControlNodePtr = S.NodeGetNodePtr(S.NodeCreateDOFs({Node::eDof::RELATIVEHUMIDITY}));

    int boundaryElementsGroup =
            S.BoundaryElementsCreate(upperSideElementsGroup, upperSideNodesGroup, boundaryControlNodePtr);

    S.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                        NuTo::Constraint::Value(*boundaryControlNodePtr, [initialRH](double rTime) -> double {
                            if (rTime == 0.) // || rTime > SIMULATIONTIME / 2.0)
                                return initialRH;
                            return 0.3;
                        }));

    std::vector<int> boundaryElementIDsVector;
    S.ElementGroupGetMembers(boundaryElementsGroup, boundaryElementIDsVector);


    for (int elementId : boundaryElementIDsVector)
    {
        S.ElementGetElementPtr(elementId)->SetIntegrationType(
                *S.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip));
    }


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int visualizeGroup = elementGroupMatrixID;
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DAMAGE);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::LOCAL_EQ_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);

    int visualizeGroup2 = elementGroupAggregatesID;
    S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);


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
    DCS("Damage");
    return 0;
}
