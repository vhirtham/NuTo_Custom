#include "base/tools/GNUPlot.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/mesh/MeshGenerator.h"

#include <fstream>
#include <functional>

using namespace NuTo;
using namespace NuTo::Constitutive;
using namespace NuTo::Interpolation;
using namespace NuTo::Node;

#define NUMELEMENTSPERDIRECTION 50

#define TIMESTEP 2.5
#define TIMESTEP_MIN TIMESTEP / 100
#define TIMESTEP_MAX TIMESTEP * 10
#define TIMESTEP_POSTPROCESS 2.5
#define SIMULATIONTIME 1000.0


#define NUMERICALTOLERANCE_DISP 1.e-3
#define NUMERICALTOLERANCE_MT 1e-6
#define MAXITERATION 20

#define EXTERNALFORCE 1.e9
#define TOTALYOUNGSMODULUS 2.e3


void RunShrinkageTest(double YoungsModulus, double poissonRatio, bool stressBased)
{

    Structure S(1);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);

    constexpr double specimenLength = 300.0;


    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //    auto meshData = S.ImportFromGmsh("ARAMIS.msh");

    auto meshData = MeshGenerator::Grid(
            S, (Eigen::VectorXd(1) << 0.).finished(), (Eigen::VectorXd(1) << specimenLength).finished(),
            (Eigen::VectorXi(1) << NUMELEMENTSPERDIRECTION - 1).finished(), Interpolation::eShapeType::TRUSS1D);
    int elementGroupID = meshData.first;
    int interpolationTypeID = meshData.second;


    // Section
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.ElementTotalSetSection(SectionTruss::Create(1.0));


    // Constitutive law
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix - mechanics
    int lawLE_ID = S.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    S.ConstitutiveLawSetParameterDouble(lawLE_ID, eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    S.ConstitutiveLawSetParameterDouble(lawLE_ID, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);


    // matrix - moisture transport

    int lawMoistureID = S.ConstitutiveLawCreate(eConstitutiveType::MOISTURE_TRANSPORT);
    S.ConstitutiveLawSetParameterBool(lawMoistureID, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,
                                      false);
    S.ConstitutiveLawSetParameterBool(lawMoistureID, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,
                                      false);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH,
                                        1.0e3);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV,
                                        1.0e10);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,
                                        1.0e-1);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV, 1.e5);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::DIFFUSION_EXPONENT_RH, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::DIFFUSION_EXPONENT_WV, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID,
                                        eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION, 0.56);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID,
                                        eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION, 0.26);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::MASS_EXCHANGE_RATE, 2.e-3);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::PORE_VOLUME_FRACTION, 0.25);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR, 0.0173);


    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMoistureID, eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,
            (Eigen::Vector4d() << 0., 0.19692057340725558, -0.28253538941816925, 0.22661481601091368).finished());
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMoistureID, eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,
            (Eigen::Vector4d() << 0., 0.26719233184420238, -0.41030868184510738, 0.32511635000090505).finished());

    // Strain based shrinkage
    int lawCSStra_ID = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    S.ConstitutiveLawSetParameterDouble(lawCSStra_ID, eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS, 30.e9);
    S.ConstitutiveLawSetParameterDouble(lawCSStra_ID, eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS, 7.5e9);
    S.ConstitutiveLawSetParameterDouble(lawCSStra_ID, eConstitutiveParameter::TEMPERATURE, 293);

    // Stress based shrinkage
    int lawCSStre_ID = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED);
    S.ConstitutiveLawSetParameterDouble(lawCSStre_ID, eConstitutiveParameter::TEMPERATURE, 293);
    int lawSAWSStre_ID = S.ConstitutiveLawCreate(eConstitutiveType::SHRINKAGE_SAW_STRESS_BASED);
    S.ConstitutiveLawSetParameterDouble(lawCSStre_ID, eConstitutiveParameter::TEMPERATURE, 293);

    // matrix

    int lawAOID = S.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);
    int lawAIE_ID = S.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);

    auto lawAOPtr = static_cast<AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAOID));
    auto lawAIE = static_cast<AdditiveInputExplicit*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_ID));


    if (stressBased)
    {
        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawLE_ID));
        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCSStre_ID));
        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawSAWSStre_ID));
        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID));
        S.ElementGroupSetConstitutiveLaw(elementGroupID, lawAOID);
    }
    else
    {
        lawAIE->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawLE_ID));
        lawAIE->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCSStra_ID), eInput::ENGINEERING_STRAIN);

        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAIE_ID));
        lawAOPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID));
        S.ElementGroupSetConstitutiveLaw(elementGroupID, lawAOID);
    }

    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationTypeID, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);


    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationTypeID, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeID, Node::eDof::RELATIVEHUMIDITY, Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeID, Node::eDof::WATERVOLUMEFRACTION,
                           Interpolation::eTypeOrder::EQUIDISTANT1);


    S.ElementTotalConvertToInterpolationType();

    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, specimenLength, specimenLength);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));
    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double initialRH = 1.0;
    int SideNodesGroup = S.GroupCreateNodeGroup();
    int SideElementsGroup = S.GroupCreateElementGroup();

    S.GroupAddNodeCoordinateRange(SideNodesGroup, eDirection::X, 0.0 - 1.e-6, 0.0 + 1.e-6);
    S.GroupAddNodeCoordinateRange(SideNodesGroup, eDirection::X, specimenLength - 1.e-6, specimenLength + 1.e-6);
    S.GroupAddElementsFromNodes(SideElementsGroup, SideNodesGroup, false);

    auto boundaryControlNodePtr = S.NodeGetNodePtr(S.NodeCreateDOFs({Node::eDof::RELATIVEHUMIDITY}));

    int boundaryElementsGroup = S.BoundaryElementsCreate(SideElementsGroup, SideNodesGroup, boundaryControlNodePtr);

    S.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                        NuTo::Constraint::Value(*boundaryControlNodePtr, [initialRH](double rTime) -> double {
                            if (rTime == 0.)
                                return initialRH;
                            return 0.4;
                        }));

    std::vector<int> boundaryElementIDsVector;
    S.ElementGroupGetMembers(boundaryElementsGroup, boundaryElementIDsVector);


    for (int elementId : boundaryElementIDsVector)
    {
        S.ElementGetElementPtr(elementId)->SetIntegrationType(
                *S.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType0DBoundary));
    }


    // Initial Values
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    int allElementsGroup_Id = S.GroupCreateElementGroup();
    S.GroupAddElementsTotal(allElementsGroup_Id);

    std::vector<int> nodeIDVectorMatrix;
    S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(allElementsGroup_Id), nodeIDVectorMatrix);

    double initialWV = S.ConstitutiveLawGetEquilibriumWaterVolumeFraction(
            lawMoistureID, initialRH,
            S.ConstitutiveLawGetParameterFullVectorDouble(
                    lawMoistureID, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
    for (auto nodeID : nodeIDVectorMatrix)
    {
        auto nodePtr = S.NodeGetNodePtr(nodeID);
        nodePtr->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, initialRH);
        nodePtr->Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, initialWV);
    }

    // Initialize Static Data
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto lawMoisturePtr = S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID);
    std::vector<int> elementIDsMatrixVector;
    S.ElementGroupGetMembers(allElementsGroup_Id, elementIDsMatrixVector);

    for (unsigned int i = 0; i < elementIDsMatrixVector.size(); i++)
    {
        for (int theIP = 0; theIP < S.ElementGetElementPtr(elementIDsMatrixVector[i])->GetNumIntegrationPoints();
             theIP++)
        {
            NuTo::Constitutive::IPAdditiveOutput& ipLawAO = dynamic_cast<NuTo::Constitutive::IPAdditiveOutput&>(
                    S.ElementGetElementPtr(i)->GetIPData().GetIPConstitutiveLaw(theIP));


            NuTo::Constitutive::StaticData::DataMoistureTransport& moistureData =
                    ipLawAO.GetSublawData<NuTo::MoistureTransport>(lawMoisturePtr).GetData(); // finally the data.

            moistureData.SetLastSorptionCoeff(lawMoisturePtr->GetParameterFullVectorDouble(
                    NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
            moistureData.SetCurrentSorptionCoeff(lawMoisturePtr->GetParameterFullVectorDouble(
                    NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
            moistureData.SetLastRelHumValue(initialRH);
            moistureData.SetDesorption(true);
        }
    }


    // Loads
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const auto& rightNode = S.NodeGetAtCoordinate((Eigen::VectorXd(1) << specimenLength).finished());


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //    int visualizeGroup = S.GroupCreate(eGroupId::Elements);
    //    S.GroupAddElementsTotal(visualizeGroup);

    int visualizeGroup = elementGroupID;

    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
    // S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::SHRINKAGE_STRAIN);


    // Set result directory
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    std::string resultDir = "1DCreepShrinkageTestResults";

    // Setup custom timestepping
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    NM.SetAutomaticTimeStepping(true);
    NM.GetTimeControl().SetTimeStepFunction(
            [](TimeControl& timeControl, int iterations, int maxIterations, bool converged) -> double {
                if (converged)
                {
                    if (iterations < 0.25 * maxIterations)
                        return timeControl.GetTimeStep() * 2.0;
                    else
                        return timeControl.GetTimeStep();
                }
                else
                {
                    timeControl.RestorePreviousTime();
                    return timeControl.GetTimeStep() * 0.5;
                }
            });

    NM.GetTimeControl().AdjustTimestep(0, 1, true);

    // Setup custom postprocessing
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Plot::GNUPlot gnuplot(true);

    const int numTimesteps = std::ceil(SIMULATIONTIME / TIMESTEP_POSTPROCESS) + 1;

    int lastCallbackTime = -1;
    double nextPostprocessTime = 0;


    std::vector<double> timeVector;
    std::vector<double> dispVector;


    ElementBase* someElement = S.ElementGetElementPtr(0);
    Eigen::MatrixXd elementStress;

    NM.PostProcessing().SetCallback([&](const StructureBase& S, const TimeControl& TC) {
        assert(lastCallbackTime < TC.GetCurrentTime());
        lastCallbackTime = TC.GetCurrentTime();


        if (TC.GetCurrentTime() >= nextPostprocessTime || TC.GetCurrentTime() == 0.)
        {
            nextPostprocessTime += TIMESTEP_POSTPROCESS;
            int index = std::floor(TC.GetCurrentTime() / TIMESTEP_POSTPROCESS);
            assert(index < numTimesteps);

            someElement->GetIntegratedStress(elementStress);
            timeVector.push_back(TC.GetCurrentTime());
            dispVector.push_back(rightNode.Get(Node::eDof::DISPLACEMENTS)[0]);
            // auto elementStress = S.ElementGetEngineeringStrain(0);

            gnuplot.Clear();

            gnuplot.AddPlot(timeVector, dispVector, {255, 0, 0}, Plot::eLineType::LINES, "numerical solution");
            gnuplot.Show();
        }
    });
    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //    S.SetNumProcessors(4);
    //    S.CalculateMaximumIndependentSets();

    NM.AddCalculationStep({NuTo::Node::eDof::RELATIVEHUMIDITY, NuTo::Node::eDof::WATERVOLUMEFRACTION});
    NM.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    NM.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, NUMERICALTOLERANCE_DISP);
    NM.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY, NUMERICALTOLERANCE_MT);
    NM.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION, NUMERICALTOLERANCE_MT);

    NM.SetMaxTimeStep(TIMESTEP_MAX);
    NM.SetMinTimeStep(TIMESTEP_MIN);
    NM.SetTimeStep(TIMESTEP);

    NM.SetPerformLineSearch(false);
    NM.SetMaxNumIterations(MAXITERATION);
    NM.PostProcessing().SetMinTimeStepPlot(TIMESTEP_POSTPROCESS);
    NM.PostProcessing().SetResultDirectory(resultDir, true);
    NM.Solve(SIMULATIONTIME);
}


int main(int argc, char* argv[])
{


    RunShrinkageTest(30.e3, 0.0, true);

    return 0;
}
