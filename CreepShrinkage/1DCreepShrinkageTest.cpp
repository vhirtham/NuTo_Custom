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
using namespace NuTo::Constraint;

#define NUMELEMENTSPERDIRECTION 3

#define INITIALTIMESTEP 1.e-9
#define TIMESTEP 2500
#define TIMESTEP_MIN TIMESTEP / 100
#define TIMESTEP_MAX TIMESTEP * 10
#define TIMESTEP_POSTPROCESS 25000.0
#define SIMULATIONTIME 50000000.0


#define NUMERICALTOLERANCE_DISP 1.e-3
#define NUMERICALTOLERANCE_MT 1e-17
#define MAXITERATION 20

#define EXTERNALFORCE 1.e9
#define TOTALYOUNGSMODULUS 2.e9


void RunCreepTestl(double YoungsModulus, Eigen::VectorXd kelvinChainStiffness,
                   Eigen::VectorXd kelvinChainRetardationTimes, double poissonRatio, double extForceDisp, bool useForce,
                   Eigen::VectorXd& timeVector, Eigen::VectorXd& strainVector, Eigen::VectorXd& stressVector)
{
    assert(kelvinChainStiffness.rows() == kelvinChainRetardationTimes.rows());
    assert(kelvinChainStiffness.cols() == kelvinChainRetardationTimes.cols());

    Structure S(1);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);


    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //    auto meshData = S.ImportFromGmsh("ARAMIS.msh");

    auto meshData =
            MeshGenerator::Grid(S, (Eigen::VectorXd(1) << 0.).finished(), (Eigen::VectorXd(1) << 0.2).finished(),
                                (Eigen::VectorXi(1) << 2).finished(), Interpolation::eShapeType::TRUSS1D);
    int elementGroupID = meshData.first;
    int interpolationTypeID = meshData.second;


    // Section
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.ElementTotalSetSection(SectionTruss::Create(1.0));


    // Constitutive law
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            lawMoistureID, Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 0.0e-6);
    S.ConstitutiveLawSetParameterDouble(
            lawMoistureID, Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 0.0e-1);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,
                                        3.9e-6);
    S.ConstitutiveLawSetParameterDouble(lawMoistureID, Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV,
                                        1.17e-0);
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
                                        30.e9);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageID, Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                        7.5e9);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageID, Constitutive::eConstitutiveParameter::TEMPERATURE, 293);

    // matrix

    int lawMatrixID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int lawAdditiveInputExplicitID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);

    auto lawMatrixPtr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawMatrixID));
    auto lawAdditiveInputExplicitPtr = static_cast<NuTo::AdditiveInputExplicit*>(
            S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));

    lawAdditiveInputExplicitPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCreepID));
    // lawAdditiveInputExplicitPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawShrinkageID),
    //                                                 NuTo::Constitutive::eInput::ENGINEERING_STRAIN);


    lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));
    lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID));
    S.ElementGroupSetConstitutiveLaw(elementGroupID, lawMatrixID);


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


    // Initial Values
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    std::vector<int> nodeIDVectorMatrix;
    S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(elementGroupID), nodeIDVectorMatrix);
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
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto lawMoisturePtr = S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID);
    std::vector<int> elementIDsMatrixVector;
    S.ElementGroupGetMembers(elementGroupID, elementIDsMatrixVector);

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


    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.2, 0.2);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(rightNodesGroup, {eDirection::X}));
    auto FNCTimeDepConst = [](double time) -> double {
        if (time > 0.0)
            return -0.04;
        else
            return 0.0;
    };
    if (!useForce)
    {
        S.Constraints().Add(Node::eDof::DISPLACEMENTS,
                            Constraint::Component(rightNodesGroup, {eDirection::X}, FNCTimeDepConst));
    }


    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int SideNodesGroup = S.GroupCreateNodeGroup();
    int SideElementsGroup = S.GroupCreateElementGroup();

    S.GroupAddNodeCoordinateRange(SideNodesGroup, eDirection::X, 0.0 - 1.e-6, 0.0 + 1.e-6);
    S.GroupAddNodeCoordinateRange(SideNodesGroup, eDirection::X, 0.2 - 1.e-6, 0.2 + 1.e-6);
    S.GroupAddElementsFromNodes(SideElementsGroup, SideNodesGroup, false);

    auto boundaryControlNodePtr = S.NodeGetNodePtr(S.NodeCreateDOFs({Node::eDof::RELATIVEHUMIDITY}));

    int boundaryElementsGroup = S.BoundaryElementsCreate(SideElementsGroup, SideNodesGroup, boundaryControlNodePtr);

    S.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                        NuTo::Constraint::Value(*boundaryControlNodePtr, [initialRH](double rTime) -> double {
                            if (rTime == 0.)
                                return initialRH;
                            return initialRH;
                        }));

    std::vector<int> boundaryElementIDsVector;
    S.ElementGroupGetMembers(boundaryElementsGroup, boundaryElementIDsVector);


    for (int elementId : boundaryElementIDsVector)
    {
        S.ElementGetElementPtr(elementId)->SetIntegrationType(
                *S.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType0DBoundary));
    }


    // Loads
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const auto& rightNode = S.NodeGetAtCoordinate((Eigen::VectorXd(1) << 0.2).finished());

    Eigen::MatrixXd timeDependentLoad(3, 2);
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = INITIALTIMESTEP;
    timeDependentLoad(2, 0) = SIMULATIONTIME;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = extForceDisp;
    timeDependentLoad(2, 1) = extForceDisp;


    Eigen::VectorXd loadDirection = Eigen::VectorXd::Zero(1);
    loadDirection[ToComponentIndex(eDirection::X)] = 1;
    if (useForce)
    {
        int load = S.LoadCreateNodeForce(&rightNode, loadDirection, 1);
        NM.SetTimeDependentLoadCase(load, timeDependentLoad);
    }
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


    timeVector = Eigen::VectorXd::Zero(numTimesteps);
    strainVector = Eigen::VectorXd::Zero(numTimesteps);
    stressVector = Eigen::VectorXd::Zero(numTimesteps);
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
            timeVector[index] = TC.GetCurrentTime();
            strainVector[index] = rightNode.Get(Node::eDof::DISPLACEMENTS)[0] / 0.2;
            stressVector[index] = elementStress(0, 0);
            // auto elementStress = S.ElementGetEngineeringStrain(0);

            gnuplot.Clear();

            gnuplot.AddPlot(timeVector.head(index), strainVector.head(index) * -1, {255, 0, 0}, Plot::eLineType::LINES,
                            "numerical solution");
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
    Eigen::VectorXd timeVector1;
    Eigen::VectorXd strainVector1;
    Eigen::VectorXd stressVector1;
    Eigen::VectorXd timeVector2;
    Eigen::VectorXd strainVector2;
    Eigen::VectorXd stressVector2;
    Eigen::VectorXd timeVector3;
    Eigen::VectorXd strainVector3;
    Eigen::VectorXd stressVector3;

    RunCreepTestl(30.e9, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(), 0.0,
                  -3.e9, false, timeVector1, strainVector1, stressVector1);

    Plot::GNUPlot gnuplot(true);
    gnuplot.AddPlot(timeVector1, stressVector1 * -1, {255, 0, 0}, Plot::eLineType::LINES, "stiff spring");
    // gnuplot.AddPlot(timeVector2, (strainVector2) * -1, {0, 255, 0}, Plot::eLineType::LINES, "creep ");
    // gnuplot.AddPlot(timeVector3, strainVector3 * -1, {0, 0, 255}, Plot::eLineType::LINES, "soft spring");
    gnuplot.Show();
    std::ofstream file;
    file.open("1DCreepDisp.csv");
    for (unsigned int i = 0; i < timeVector1.rows(); ++i)
        file << timeVector1[i] << ", " << stressVector1[i] << ", " << strainVector1[i] << "\n";
    file.close();

    return 0;
}
