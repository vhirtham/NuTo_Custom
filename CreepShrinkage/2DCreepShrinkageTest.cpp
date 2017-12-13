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
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "visualize/VisualizeEnum.h"

#include <functional>

using namespace NuTo;
using namespace NuTo::Constraint;

#define NUMELEMENTSPERDIRECTION 3

#define INITIALTIMESTEP 1.e-9
#define TIMESTEP 25000
#define TIMESTEP_MIN TIMESTEP / 100
#define TIMESTEP_MAX TIMESTEP * 10
#define TIMESTEP_POSTPROCESS 250000.0
#define SIMULATIONTIME 25000000.0


#define NUMERICALTOLERANCE_DISP 1.e-3
#define NUMERICALTOLERANCE_MT 1e-17
#define MAXITERATION 20

#define EXTERNALFORCE 0
#define TOTALYOUNGSMODULUS 2.e9


void RunCreepTestl(double YoungsModulus, Eigen::VectorXd kelvinChainStiffness,
                   Eigen::VectorXd kelvinChainRetardationTimes, double poissonRatio, double extLoad, bool IsHomogeneous,
                   bool useShrinkage, Eigen::VectorXd& timeVector, Eigen::VectorXd& strainVector,
                   std::string simulationName)
{
    assert(kelvinChainStiffness.rows() == kelvinChainRetardationTimes.rows());
    assert(kelvinChainStiffness.cols() == kelvinChainRetardationTimes.cols());

    Structure S(2);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);


    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto meshData = S.ImportFromGmsh("Meso.msh");

    int elementGroupMatrixID = meshData[0].first;
    int interpolationTypeMatrixID = meshData[0].second;

    int elementGroupAggregatesID = meshData[1].first;
    int interpolationTypeAggregatesID = meshData[1].second;

    int allElementsGroup = S.GroupCreateElementGroup();
    S.GroupAddElementsTotal(allElementsGroup);


    // Section
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.ElementTotalSetSection(SectionPlane::Create(1.0, false));


    // Constitutive law
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix - mechanics
    int lawCreepID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::CREEP);

    S.ConstitutiveLawSetParameterDouble(lawCreepID, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                        YoungsModulus);
    S.ConstitutiveLawSetParameterDouble(lawCreepID, Constitutive::eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawCreepID, Constitutive::eConstitutiveParameter::KELVIN_CHAIN_STIFFNESS, kelvinChainStiffness);
    S.ConstitutiveLawSetParameterFullVectorDouble(lawCreepID,
                                                  Constitutive::eConstitutiveParameter::KELVIN_CHAIN_RETARDATIONTIME,
                                                  kelvinChainRetardationTimes);


    // matrix - moisture transport

    int lawMoistureID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    S.ConstitutiveLawSetParameterBool(
            lawMoistureID, Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterBool(
            lawMoistureID, Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);

    S.ConstitutiveLawSetParameterDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 1.0e-9);
    S.ConstitutiveLawSetParameterDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 1.0e-5);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,
                                        3.9e-12);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV,
                                        1.17e-7);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH,
                                        1.0);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV,
                                        1.0);
    S.ConstitutiveLawSetParameterDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION, 0.56);
    S.ConstitutiveLawSetParameterDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION, 0.26);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE, 0.e-8);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION,
                                        0.25);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID,
                                        Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR, 0.0173);


    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,
            (Eigen::Vector4d() << 0., 0.19692057340725558, -0.28253538941816925, 0.22661481601091368).finished());
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,
            (Eigen::Vector4d() << 0., 0.26719233184420238, -0.41030868184510738, 0.32511635000090505).finished());

    // Strain based shrinkage
    int lawShrinkageID =
            S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageID, Constitutive::eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,
                                        YoungsModulus * 4. / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageID, Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                        YoungsModulus * 2. / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageID, Constitutive::eConstitutiveParameter::TEMPERATURE, 293);

    // matrix

    int lawMatrixID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int lawAdditiveInputExplicitID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);

    auto lawMatrixPtr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawMatrixID));
    auto lawAdditiveInputExplicitPtr = static_cast<NuTo::AdditiveInputExplicit*>(
            S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));

    lawAdditiveInputExplicitPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCreepID));
    lawAdditiveInputExplicitPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawShrinkageID),
                                                    NuTo::Constitutive::eInput::ENGINEERING_STRAIN);


    if (useShrinkage)
        lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));
    else
        lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCreepID));
    lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID));
    S.ElementGroupSetConstitutiveLaw(elementGroupMatrixID, lawMatrixID);


    // aggregates

    int lawLinearElasticID =
            S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(lawLinearElasticID, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                        YoungsModulus * 4);
    S.ConstitutiveLawSetParameterDouble(lawLinearElasticID, Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                        poissonRatio);

    if (!IsHomogeneous)
        S.ElementGroupSetConstitutiveLaw(elementGroupAggregatesID, lawLinearElasticID);
    else
        S.ElementGroupSetConstitutiveLaw(elementGroupAggregatesID, lawMatrixID);


    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationTypeMatrixID,
                                          NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);


    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::DISPLACEMENTS,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::RELATIVEHUMIDITY,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::WATERVOLUMEFRACTION,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    // aggregates
    if (!IsHomogeneous)
        S.InterpolationTypeAdd(interpolationTypeAggregatesID, Node::eDof::DISPLACEMENTS,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
    else
    {
        S.InterpolationTypeAdd(interpolationTypeAggregatesID, Node::eDof::DISPLACEMENTS,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
        S.InterpolationTypeAdd(interpolationTypeAggregatesID, Node::eDof::RELATIVEHUMIDITY,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
        S.InterpolationTypeAdd(interpolationTypeAggregatesID, Node::eDof::WATERVOLUMEFRACTION,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
    }

    S.ElementTotalConvertToInterpolationType();


    // Initial Values
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    std::vector<int> nodeIDVectorMatrix;
    if (IsHomogeneous)
    {
        S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(allElementsGroup), nodeIDVectorMatrix);
    }
    else
        S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(elementGroupMatrixID), nodeIDVectorMatrix);
    double initialRH = 0.95;
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
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto lawMoisturePtr = S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID);
    std::vector<int> elementIDsMatrixVector;
    if (IsHomogeneous)
        S.ElementGroupGetMembers(allElementsGroup, elementIDsMatrixVector);
    else
        S.ElementGroupGetMembers(elementGroupMatrixID, elementIDsMatrixVector);

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


    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto virtualNodePtr = S.NodeGetNodePtr(S.NodeCreate(Eigen::VectorXd::Ones(2) * -1, {Node::eDof::DISPLACEMENTS}));

    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.2, 0.2);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));
    for (auto& itNode : rightNodesGroup)
    {
        S.Constraints().Add(Node::eDof::DISPLACEMENTS,
                            Equation({Term(*virtualNodePtr, ToComponentIndex(eDirection::X), 1),
                                      Term(*S.NodeGetNodePtr(itNode.first), ToComponentIndex(eDirection::X), -1)}));
    }

    // Additional 2D constraints
    auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(2));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(*virtualNodePtr, {eDirection::Y}));


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


    // Loads
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigen::MatrixXd timeDependentLoad(3, 2);
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = INITIALTIMESTEP;
    timeDependentLoad(2, 0) = SIMULATIONTIME;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = extLoad;
    timeDependentLoad(2, 1) = extLoad;


    Eigen::VectorXd loadDirection = Eigen::VectorXd::Zero(2);
    loadDirection[ToComponentIndex(eDirection::X)] = 1;
    int load = S.LoadCreateNodeForce(virtualNodePtr, loadDirection, 1);
    NM.SetTimeDependentLoadCase(load, timeDependentLoad);


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //    int visualizeGroup = S.GroupCreate(eGroupId::Elements);
    //    S.GroupAddElementsTotal(visualizeGroup);


    int visualizeGroup = elementGroupMatrixID;
    if (IsHomogeneous)
        visualizeGroup = allElementsGroup;

    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
    if (useShrinkage)
        S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::SHRINKAGE_STRAIN);


    if (!IsHomogeneous)
    {
        int visualizeGroup2 = elementGroupAggregatesID;
        S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::DISPLACEMENTS);
        S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::ENGINEERING_STRAIN);
        S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::ENGINEERING_STRESS);
        S.AddVisualizationComponent(visualizeGroup2, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    }


    // Setup custom timestepping
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

    Plot::GNUPlot gnuplot;

    const int numTimesteps = std::ceil(SIMULATIONTIME / TIMESTEP_POSTPROCESS) + 1;

    int lastCallbackTime = -1;
    double nextPostprocessTime = 0;

    const auto& rightLowerNode = S.NodeGetAtCoordinate((Eigen::VectorXd(2) << 0.2, 0.2).finished());
    timeVector = Eigen::VectorXd::Zero(numTimesteps);
    strainVector = Eigen::VectorXd::Zero(numTimesteps);

    NM.PostProcessing().SetCallback([&](const StructureBase& S, const TimeControl& TC) {
        assert(lastCallbackTime < TC.GetCurrentTime());
        lastCallbackTime = TC.GetCurrentTime();


        if (TC.GetCurrentTime() >= nextPostprocessTime || TC.GetCurrentTime() == 0.)
        {
            nextPostprocessTime += TIMESTEP_POSTPROCESS;
            int index = std::floor(TC.GetCurrentTime() / TIMESTEP_POSTPROCESS);
            assert(index < numTimesteps);
            timeVector[index] = TC.GetCurrentTime();
            strainVector[index] = rightLowerNode.Get(Node::eDof::DISPLACEMENTS)[0] / 0.2;

            gnuplot.Clear();

            gnuplot.AddPlot(timeVector.head(index), strainVector.head(index) * -1, {255, 0, 0}, Plot::eLineType::LINES,
                            "numerical solution");
            gnuplot.Show();
        }
    });

    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.SetNumProcessors(4);
    S.CalculateMaximumIndependentSets();

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
    NM.PostProcessing().SetResultDirectory(simulationName, true);
    NM.Solve(SIMULATIONTIME);
}


void old()
{
    Structure S(2);
    auto mshData = S.ImportFromGmsh("ARAMIS.msh");


    // create section
    double thickness = 20.0;
    auto section = SectionPlane::Create(thickness, false);
    S.ElementTotalSetSection(section);

    // group IDs
    int GIMatrix = mshData[0].first;
    int GIAggregate = mshData[1].first;

    // interpolation type IDs
    int ITIMatrix = mshData[0].second;
    int ITIAggregate = mshData[1].second;

    // Matrix
    int CILEMatrix = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(CILEMatrix, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30 * 10e9);
    S.ConstitutiveLawSetParameterDouble(CILEMatrix, Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);
    S.ElementGroupSetConstitutiveLaw(GIMatrix, CILEMatrix);
    S.InterpolationTypeAdd(ITIMatrix, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);

    // Aggregates
    int CILEAggregates = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(CILEAggregates, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                        60 * 10e9);
    S.ConstitutiveLawSetParameterDouble(CILEAggregates, Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);
    S.ElementGroupSetConstitutiveLaw(GIAggregate, CILEAggregates);
    S.InterpolationTypeAdd(ITIAggregate, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);

    S.ElementTotalConvertToInterpolationType();


    // Groups
    auto GLowBound = S.GroupGetNodeCoordinateRange(eDirection::Y, 0.0, 0.0);
    auto GUpBound = S.GroupGetNodeCoordinateRange(eDirection::Y, 0.2, 0.2);
    auto GILowLeft = S.GroupCreateNodeGroup();
    const auto& GLowLeft = *(S.GroupGetGroupPtr(GILowLeft)->AsGroupNode());
    S.GroupAddNodeFunction(GILowLeft, [](NodeBase* node) -> bool {
        if (node->Get(NuTo::Node::eDof::COORDINATES)[0] == 0.0 && node->Get(NuTo::Node::eDof::COORDINATES)[1] == 0.0)
            return true;
        else
            return false;
    });

    assert(GLowBound.GetNumMembers() > 0);
    assert(GUpBound.GetNumMembers() > 0);
    assert(GLowLeft.GetNumMembers() == 1);


    // Constraints
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(GLowLeft, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(GLowBound, {eDirection::Y}));


    auto FNCTimeDepConst = [](double time) -> double {
        if (time > 0.0)
            return 1.0;
        else
            return 0.0;
    };

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(GUpBound, {eDirection::Y}, FNCTimeDepConst));

    // Loads

    //    S.LoadCreateNodeGroupForce(GUpBound);


    // Visualization
    int GIVis = S.GroupCreate(eGroupId::Elements);
    S.GroupAddElementsTotal(GIVis);

    S.AddVisualizationComponent(GIVis, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(GIVis, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(GIVis, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(GIVis, eVisualizeWhat::SHRINKAGE_STRAIN);
    S.AddVisualizationComponent(GIVis, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);


    // Solve

    NewmarkDirect NM(&S);
    NM.SetAutomaticTimeStepping(false);
    NM.SetTimeStep(1800.0);
    NM.SetPerformLineSearch(false);
    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-2);

    NM.PostProcessing().SetResultDirectory("ARAMISResults", true);
    NM.Solve(3600.0);
}


int main(int argc, char* argv[])
{

    Eigen::VectorXd timeVector1;
    Eigen::VectorXd strainVector1;
    Eigen::VectorXd timeVector2;
    Eigen::VectorXd strainVector2;
    Eigen::VectorXd timeVector3;
    Eigen::VectorXd strainVector3;

    //    RunCreepTestl(3.e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, false, "CSS_elastic_meso");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 1000000.).finished(),
    //    0.2,
    //                  false, "CSS_soft_meso");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(),
    //    0.2,
    //                  false, "CSS_stiff_meso");

    //    RunCreepTestl(3.e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, true,
    //    "CSS_elastic_homogeneous");
    RunCreepTestl(30.e9, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(), 0.2,
                  -3e9 * 0.2, false, false, timeVector1, strainVector1, "CSS_meso_noshrinkage");
    //    RunCreepTestl(30.e9, (Eigen::VectorXd::Zero(0)), (Eigen::VectorXd::Zero(0)), 0.2, false, "CSS_meso_nocreep");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(),
    //    0.2,
    //                  true, "CSS_stiff_homogeneous");

    Plot::GNUPlot gnuplot(true);
    gnuplot.AddPlot(timeVector1, strainVector1 * -1, {255, 0, 0}, Plot::eLineType::LINES, "2D Creep Meso");
    // gnuplot.AddPlot(timeVector2, (strainVector2) * -1, {0, 255, 0}, Plot::eLineType::LINES, "creep ");
    // gnuplot.AddPlot(timeVector3, strainVector3 * -1, {0, 0, 255}, Plot::eLineType::LINES, "soft spring");
    gnuplot.Show();
    std::ofstream file;
    file.open("2DCreepMeso.csv");
    for (unsigned int i = 0; i < timeVector1.rows(); ++i)
        file << timeVector1[i] << ", " << strainVector1[i] << "\n";
    file.close();
    return 0;
}
