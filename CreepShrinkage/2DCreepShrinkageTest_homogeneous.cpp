#include "base/tools/GNUPlot.h"
#include "base/tools/PlotEnum.h"
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
#include "mechanics/mesh/MeshGenerator.h"
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
                   Eigen::VectorXd kelvinChainRetardationTimes, double poissonRatio, double extForce,
                   bool useMoistureTransport, bool useStrainBasedShrinkage, bool useCreepConstraints,
                   Eigen::VectorXd& timeVector, Eigen::VectorXd& strainVector,
                   std::string simulationName = "Creep_Shrinkage_Simulation")
{
    assert(kelvinChainStiffness.rows() == kelvinChainRetardationTimes.rows());
    assert(kelvinChainStiffness.cols() == kelvinChainRetardationTimes.cols());

    Structure S(2);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);


    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto meshData = MeshGenerator::Grid(S, (Eigen::VectorXd(2) << 0., 0.).finished(),
                                        (Eigen::VectorXd(2) << 0.2, 0.2).finished(),
                                        (Eigen::VectorXi(2) << 20, 20).finished(), Interpolation::eShapeType::QUAD2D);


    int elementGroupMatrixID = meshData.first;
    int interpolationTypeMatrixID = meshData.second;


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

    int lawShrinkageStrainBasedID =
            S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageStrainBasedID,
                                        Constitutive::eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,
                                        YoungsModulus * 4. / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageStrainBasedID,
                                        Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                        YoungsModulus * 2. / 3.);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageStrainBasedID, Constitutive::eConstitutiveParameter::TEMPERATURE,
                                        293);

    // Stress based shrinkage

    int lawShrinkageStressBasedID =
            S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED);
    S.ConstitutiveLawSetParameterDouble(lawShrinkageStressBasedID, Constitutive::eConstitutiveParameter::TEMPERATURE,
                                        293);

    // matrix

    int lawMatrixID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int lawAdditiveInputExplicitID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);

    auto lawMatrixPtr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(lawMatrixID));
    auto lawAdditiveInputExplicitPtr = static_cast<NuTo::AdditiveInputExplicit*>(
            S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));


    if (useStrainBasedShrinkage)
    {
        lawAdditiveInputExplicitPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCreepID));
        if (useMoistureTransport)
            lawAdditiveInputExplicitPtr->AddConstitutiveLaw(
                    *S.ConstitutiveLawGetConstitutiveLawPtr(lawShrinkageStrainBasedID),
                    NuTo::Constitutive::eInput::ENGINEERING_STRAIN);
        lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawAdditiveInputExplicitID));
    }
    else
    {
        lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawCreepID));
        if (useMoistureTransport)
            lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawShrinkageStressBasedID));
    }
    if (useMoistureTransport)
        lawMatrixPtr->AddConstitutiveLaw(*S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID));
    S.ElementGroupSetConstitutiveLaw(elementGroupMatrixID, lawMatrixID);


    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationTypeMatrixID,
                                          NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);


    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::DISPLACEMENTS,
                           Interpolation::eTypeOrder::EQUIDISTANT1);
    if (useMoistureTransport)
    {
        S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::RELATIVEHUMIDITY,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
        S.InterpolationTypeAdd(interpolationTypeMatrixID, Node::eDof::WATERVOLUMEFRACTION,
                               Interpolation::eTypeOrder::EQUIDISTANT1);
    }

    S.ElementTotalConvertToInterpolationType();


    // Initial Values
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double initialRH = 0.0;
    double initialWV = 0.0;

    // nodes
    if (useMoistureTransport)
    {
        std::vector<int> nodeIDVectorMatrix;
        S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(elementGroupMatrixID), nodeIDVectorMatrix);
        initialRH = 0.95;
        initialWV = S.ConstitutiveLawGetEquilibriumWaterVolumeFraction(
                lawMoistureID, initialRH,
                S.ConstitutiveLawGetParameterFullVectorDouble(
                        lawMoistureID, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
        for (auto nodeID : nodeIDVectorMatrix)
        {
            auto nodePtr = S.NodeGetNodePtr(nodeID);
            nodePtr->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, initialRH);
            nodePtr->Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, initialWV);
        }
    }

    // Initialize Static Data
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (useMoistureTransport)
    {
        auto lawMoisturePtr = S.ConstitutiveLawGetConstitutiveLawPtr(lawMoistureID);
        std::vector<int> elementIDsMatrixVector;

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
    }

    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    NodeBase* virtualNodePtr = nullptr;

    auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(2));
    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.2, 0.2);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    if (useCreepConstraints)
    {
        virtualNodePtr = S.NodeGetNodePtr(S.NodeCreate(Eigen::VectorXd::Ones(2) * -1, {Node::eDof::DISPLACEMENTS}));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));
        for (auto& itNode : rightNodesGroup)
        {
            S.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Equation({Term(*virtualNodePtr, ToComponentIndex(eDirection::X), 1),
                                          Term(*S.NodeGetNodePtr(itNode.first), ToComponentIndex(eDirection::X), -1)}));
        }

        // Additional 2D constraints
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(*virtualNodePtr, {eDirection::Y}));
    }
    else
    {


        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::X}));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS,
                            Constraint::Component(S.NodeGetAtCoordinate((Eigen::VectorXd(2) << 0.2, 0.0).finished()),
                                                  {eDirection::Y}));
    }

    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (useMoistureTransport)
    {
        int upperSideNodesGroup = S.GroupCreateNodeGroup();
        int upperSideElementsGroup = S.GroupCreateElementGroup();

        S.GroupAddNodeCoordinateRange(upperSideNodesGroup, eDirection::Y, 0. - 1.e-6, 0. + 1.e-6);
        S.GroupAddElementsFromNodes(upperSideElementsGroup, upperSideNodesGroup, false);

        auto boundaryControlNodePtr = S.NodeGetNodePtr(S.NodeCreateDOFs({Node::eDof::RELATIVEHUMIDITY}));

        int boundaryElementsGroup =
                S.BoundaryElementsCreate(upperSideElementsGroup, upperSideNodesGroup, boundaryControlNodePtr);

        S.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                            NuTo::Constraint::Value(*boundaryControlNodePtr, [initialRH](double rTime) -> double {
                                if (rTime == 0.) //|| rTime > SIMULATIONTIME / 2.0)
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
    }
    // Loads
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigen::MatrixXd timeDependentLoad(5, 2);
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = INITIALTIMESTEP;
    timeDependentLoad(2, 0) = SIMULATIONTIME / 2;
    timeDependentLoad(3, 0) = SIMULATIONTIME / 2 + TIMESTEP;
    timeDependentLoad(4, 0) = SIMULATIONTIME;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = extForce;
    timeDependentLoad(2, 1) = extForce;
    timeDependentLoad(3, 1) = extForce;
    timeDependentLoad(4, 1) = extForce;


    Eigen::VectorXd loadDirection = Eigen::VectorXd::Zero(2);
    loadDirection[ToComponentIndex(eDirection::X)] = 1;
    if (useCreepConstraints)
    {
        int load = S.LoadCreateNodeForce(virtualNodePtr, loadDirection, 1);
        NM.SetTimeDependentLoadCase(load, timeDependentLoad);
    }

    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    int visualizeGroup = elementGroupMatrixID;

    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    if (useMoistureTransport)
    {
        S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::RELATIVE_HUMIDITY);
        S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::WATER_VOLUME_FRACTION);
        if (useStrainBasedShrinkage)
            S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::SHRINKAGE_STRAIN);
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
            strainVector[index] = rightLowerNode.Get(Node::eDof::DISPLACEMENTS)[0];

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

    if (useMoistureTransport)
    {
        NM.AddCalculationStep({NuTo::Node::eDof::RELATIVEHUMIDITY, NuTo::Node::eDof::WATERVOLUMEFRACTION});
        NM.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY, NUMERICALTOLERANCE_MT);
        NM.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION, NUMERICALTOLERANCE_MT);
    }
    NM.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    NM.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, NUMERICALTOLERANCE_DISP);

    NM.SetMaxTimeStep(TIMESTEP_MAX);
    NM.SetMinTimeStep(TIMESTEP_MIN);
    NM.SetTimeStep(TIMESTEP);

    NM.SetPerformLineSearch(false);
    NM.SetMaxNumIterations(MAXITERATION);
    NM.PostProcessing().SetMinTimeStepPlot(TIMESTEP_POSTPROCESS);
    NM.PostProcessing().SetResultDirectory(simulationName, true);
    NM.Solve(SIMULATIONTIME);
}


int main(int argc, char* argv[])
{

    //    RunCreepTestl(3.e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, "CSS_elastic_homogeneous");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 1000000.).finished(),
    //    0.2,
    //                  "CSS_soft_homogeneous");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(),
    //    0.2,
    //                  "CSS_stiff_homogeneous");

    Eigen::VectorXd timeVector1;
    Eigen::VectorXd strainVector1;
    Eigen::VectorXd timeVector2;
    Eigen::VectorXd strainVector2;
    Eigen::VectorXd timeVector3;
    Eigen::VectorXd strainVector3;

    //    RunCreepTestl(3.e10, (Eigen::VectorXd(2) << 30.e9, 30.e6).finished(),
    //                  (Eigen::VectorXd(2) << 1000000., 100000000000.).finished(), 0.2, -15.e6, false, false,
    //                  timeVector1,
    //                  strainVector1, "CSS_stiff_homogeneous");

    //    RunCreepTestl(3.e10, (Eigen::VectorXd(2) << 30.e9, 30.e6).finished(),
    //                  (Eigen::VectorXd(2) << 1000000., 100000000000.).finished(), 0.2, 0, true, false, timeVector2,
    //                  strainVector2, "CSS_stiff_homogeneous");
    //    RunCreepTestl(3.e10, (Eigen::VectorXd(2) << 30.e9, 30.e6).finished(),
    //                  (Eigen::VectorXd(2) << 1000000., 100000000000.).finished(), 0.2, -15.e6, true, false,
    //                  timeVector3,
    //                  strainVector3, "CSS_stiff_homogeneous");

    //    Plot::GNUPlot gnuplot(true);
    //    gnuplot.AddPlot(timeVector1, strainVector1 * -1, {255, 0, 0}, Plot::eLineType::LINES, "Creep");
    //    gnuplot.AddPlot(timeVector2, (strainVector2 + strainVector1) * -1, {0, 255, 0}, Plot::eLineType::LINES,
    //                    "Creep + Shrinkage");

    //    gnuplot.AddPlot(timeVector3, strainVector3 * -1, {0, 0, 255}, Plot::eLineType::LINES,
    //                    "combined Creep and Shrinkage");
    //    gnuplot.Show();


    RunCreepTestl(3.e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, 0, true, false, true, timeVector1,
                  strainVector1, "CSS_stiff_spring");
    RunCreepTestl(3.e10, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(), 0.2, 0,
                  true, false, true, timeVector2, strainVector2, "CSS_creep");
    RunCreepTestl(1.5e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, 0, true, false, true, timeVector3,
                  strainVector3, "CSS_soft_spring");


    //    RunCreepTestl(30.e9, (Eigen::VectorXd(1) << 30.e9).finished(), (Eigen::VectorXd(1) << 10000000.).finished(),
    //    0.0,
    //                  -3e9 * 0.1, false, true, true, timeVector1, strainVector1, "test");
    //    RunCreepTestl(3.e10, Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), 0.2, 0, true, false, true,
    //    timeVector2,
    //                  strainVector2, "stress based");

    Plot::GNUPlot gnuplot(true);
    gnuplot.AddPlot(timeVector1, strainVector1 * -1, {255, 0, 0}, Plot::eLineType::LINES, "stiff");
    gnuplot.AddPlot(timeVector2, strainVector2 * -1, {0, 255, 0}, Plot::eLineType::LINES, "creep ");
    gnuplot.AddPlot(timeVector3, strainVector3 * -1, {0, 0, 255}, Plot::eLineType::LINES, "soft");
    gnuplot.Show();

    std::ofstream file;
    file.open("2DShrinkageCreep.csv");
    for (unsigned int i = 0; i < timeVector1.rows(); ++i)
        file << timeVector1[i] << ", " << strainVector1[i] << ", " << strainVector2[i] << ", " << strainVector3[i]
             << "\n";
    file.close();

    return 0;
}
