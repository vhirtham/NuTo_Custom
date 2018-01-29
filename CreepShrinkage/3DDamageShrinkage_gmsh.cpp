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
#include <iostream>
#include <fstream>

using namespace NuTo;
using namespace NuTo::Constitutive;
using namespace NuTo::Interpolation;
using namespace NuTo::Node;


// Constants
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// time
constexpr double simulationTime = 2.;
constexpr int numTimesteps = 16;
constexpr int numWriteSteps = 16;

// mechanics
constexpr double youngsModulus = 30000;
constexpr double poissonRatio = 0.2;
constexpr double tensileStrength = 4.;

// mesh
constexpr double cylinderDiameter = 100.;
constexpr double cylinderHeight = 300.;

// iteration
constexpr double tolDisp = 1e-6;
constexpr double tolNLES = 1e-6;
constexpr double tolWF = 1e-3;
constexpr double tolRH = 1e-3;
constexpr int maxIter = 10;

// Helper functions and classes
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct SimulationSetup
{
    // Drying
    double RH_env = 1.0;
    double RH_initial = 1.0;
    double dryingTime = 0.;
    double dryingTimeRamp = 40;
    double initialTimestepDrying = 1.0;
    int dryingMaxIter = 10;

    // Compression Test
    double maxCompression = 2;
    double compressionTime = 1.;
    double maxCompressionDeltaTime = 0.001;

    // Post Processing
    int numWritesDrying = 4;
    int numWritesCompression = 100;
};


double ArgToDouble(std::string Arg)
{
    try
    {
        double number = std::stod(Arg);
        return number;
    }
    catch (std::invalid_argument e)
    {
        throw Exception(__PRETTY_FUNCTION__, "Could not cast \"" + Arg + "\" into a double!");
    }
}

// Simulation
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void DS(SimulationSetup setup, bool useStressbased, std::string simulationName)
{

    // Check Setup
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (setup.dryingTime < setup.dryingTimeRamp)
        throw Exception(__PRETTY_FUNCTION__,
                        "Drying time is smaller than the time to reduce the relative humidity to its target value");
    setup.maxCompression = std::abs(setup.maxCompression);

    // Structure and TimeIntegration
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Structure S(3);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);
    S.SetVerboseLevel(0);

    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto meshData = S.ImportFromGmsh("Cylinder.msh");

    int elementGroup_Id = meshData[0].first;
    int interpolationType_Id = meshData[0].second;

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
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    S.ConstitutiveLawSetParameterDouble(lawGD_Id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, tensileStrength * 10);
    S.ConstitutiveLawSetDamageLaw(
            lawGD_Id, DamageLawExponential::Create(tensileStrength / youngsModulus, tensileStrength / 0.01));

    // moisture transport

    // matrix - moisture transport
    int lawMT_Id = S.ConstitutiveLawCreate(eConstitutiveType::MOISTURE_TRANSPORT);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterBool(lawMT_Id, eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, false);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 1.0e4);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 5.0e9);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH, 5.0e-1);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV, 1.0e5);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_RH, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_WV, 5.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION,
                                        0.56);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION,
                                        0.26);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::MASS_EXCHANGE_RATE, 1.e-2);
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

    S.InterpolationTypeSetIntegrationType(interpolationType_Id, eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip);

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
    auto& topNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::Z, 300 - 1e-6, 300 + 1e+6);
    int topElementsGroup_Id = S.GroupCreateElementGroup();
    int bottomElementsGroup_Id = S.GroupCreateElementGroup();
    S.GroupAddElementsFromNodes(topElementsGroup_Id, S.GroupGetId(&topNodesGroup), false);
    S.GroupAddElementsFromNodes(bottomElementsGroup_Id, S.GroupGetId(&bottomNodesGroup), false);


    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int BoundaryNodesGroup = S.GroupCreateNodeGroup();
    S.GroupAddNodeFunction(BoundaryNodesGroup, [](NodeBase* node) -> bool {
        auto coords = node->Get(Node::eDof::COORDINATES);
        auto r = std::sqrt(coords[0] * coords[0] + coords[1] * coords[1]);
        if (r > std::abs(cylinderDiameter / 2.) - 1e-6 || coords[2] >= cylinderHeight - 1e-3 || coords[2] <= 0 + 1e-3)
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
                        NuTo::Constraint::Value(*boundaryControlNodePtr, [&setup](double rTime) -> double {
                            if (rTime == 0.) // || rTime > SIMULATIONTIME / 2.0)
                                return setup.RH_initial;
                            else if (rTime < setup.dryingTimeRamp)
                                return setup.RH_initial -
                                       rTime / setup.dryingTimeRamp * (setup.RH_initial - setup.RH_env);
                            else
                                return setup.RH_env;
                        }));

    std::vector<int> boundaryElementIDsVector;
    S.ElementGroupGetMembers(boundaryElementsGroup, boundaryElementIDsVector);


    for (int elementId : boundaryElementIDsVector)
    {
        S.ElementGetElementPtr(elementId)->SetIntegrationType(
                *S.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType2D3NGauss16Ip));
    }

    // Initial Values
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // nodes

    std::vector<int> nodeIDVectorMatrix;

    int allElementsGroup_Id = S.GroupCreateElementGroup();
    S.GroupAddElementsTotal(allElementsGroup_Id);

    S.NodeGroupGetMembers(S.GroupCreateNodeGroupFromElements(allElementsGroup_Id), nodeIDVectorMatrix);
    double initialWV = S.ConstitutiveLawGetEquilibriumWaterVolumeFraction(
            lawMT_Id, setup.RH_initial,
            S.ConstitutiveLawGetParameterFullVectorDouble(
                    lawMT_Id, Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
    for (auto nodeID : nodeIDVectorMatrix)
    {
        auto nodePtr = S.NodeGetNodePtr(nodeID);
        if (nodePtr->IsDof(NuTo::Node::eDof::RELATIVEHUMIDITY))
            nodePtr->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, setup.RH_initial);
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
            moistureData.SetLastRelHumValue(setup.RH_initial);
            moistureData.SetDesorption(true);
        }
    }
    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    assert(bottomNodesGroup.GetNumMembers() > 0);
    assert(topNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Z}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::Z},
                                                                         [&setup](double t) {
                                                                             if (t < setup.dryingTime)
                                                                                 return 0.0;
                                                                             return (t - setup.dryingTime) /
                                                                                    (setup.compressionTime) *
                                                                                    (-1 * setup.maxCompression);
                                                                         }));
    auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(3), 1e-6);
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::X}));
    //    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::Y}));

    // Postprocessing for compression test
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    double topElementsVolume = S.ElementGroupGetVolume(topElementsGroup_Id);
    double bottomElementsVolume = S.ElementGroupGetVolume(bottomElementsGroup_Id);
    //    S.ElementGroupGetAverageStress(topElementsGroup_Id, S.ElementGroupGetVolume(topElementsGroup_Id),
    //    StressVector);


    std::ofstream outputFile;
    outputFile.open("SD.txt");
    outputFile << "Time, displacement_axialTop, displacement_axialBottom, sigma_axialTop, sigma_axialBottom\n";


    NM.PostProcessing().SetCallback([&S, &outputFile, &topNodesGroup, topElementsVolume, topElementsGroup_Id,
                                     &bottomNodesGroup, bottomElementsVolume,
                                     bottomElementsGroup_Id](const StructureBase& S, const TimeControl& TC) {
        outputFile << TC.GetCurrentTime() << ", ";
        double dispSum = 0;
        for (auto& iterator : topNodesGroup)
        {
            NodeBase& node = *iterator.second;
            dispSum += node.Get(Node::eDof::DISPLACEMENTS)[2];
        }
        outputFile << dispSum / topNodesGroup.size() << ", ";

        dispSum = 0;
        for (auto& iterator : bottomNodesGroup)
        {
            NodeBase& node = *iterator.second;
            dispSum += node.Get(Node::eDof::DISPLACEMENTS)[2];
        }
        outputFile << dispSum / bottomNodesGroup.size() << ", ";

        Eigen::MatrixXd StressVector = Eigen::MatrixXd::Zero(6, 1);
        S.ElementGroupGetAverageStress(topElementsGroup_Id, topElementsVolume, StressVector);
        outputFile << StressVector(2, 0) << ", ";

        Eigen::MatrixXd StressVector2 = Eigen::MatrixXd::Zero(6, 1);
        ;
        S.ElementGroupGetAverageStress(bottomElementsGroup_Id, bottomElementsVolume, StressVector2);
        outputFile << StressVector2(2, 0) << "\n";
        outputFile.flush();

    });

    // Setup custom timestepping
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    NM.GetTimeControl().SetTimeStepFunction([&setup](TimeControl& timeControl, int iterations, int maxIterations,
                                                     bool converged) -> double {
        struct TimestepHelper
        {
            double prevTimestep = 0.;
            bool isTimestepAdjusted = false;
        };

        static TimestepHelper timestepHelper;

        if (converged)
        {
            double timestep = timeControl.GetTimeStep();
            double currTime = timeControl.GetCurrentTime();

            if (timestepHelper.isTimestepAdjusted)
            {
                timestep = timestepHelper.prevTimestep;
                timestepHelper.isTimestepAdjusted = false;
            }
            else
            {
                if (iterations < 0.25 * maxIterations)
                {
                    timestep = timeControl.GetTimeStep() * 2.0;
                    if (timestep > timeControl.GetMaxTimeStep())
                        timestep = timeControl.GetMaxTimeStep();
                }
                else
                    timestep = timeControl.GetTimeStep();
            }

            if (currTime < setup.dryingTime)
            {
                double nextTime = currTime + timestep;

                double currNumWrites = std::floor(currTime / (setup.dryingTime / setup.numWritesDrying));
                double nextWriteTime = (currNumWrites + 1) * setup.dryingTime / setup.numWritesDrying;

                if ((nextTime + timeControl.GetMinTimeStep() > nextWriteTime))
                {
                    timestepHelper.prevTimestep = timestep;
                    timestepHelper.isTimestepAdjusted = true;
                    timestep = nextWriteTime - currTime +
                               1e-10; // +1e-10 to avoid floating point problems (minimal timestep)
                }
            }
            else
            {
                double nextTime = currTime + timestep;
                double currCompressionTime = currTime - setup.dryingTime;
                double currNumWrites =
                        std::floor(currCompressionTime / (setup.compressionTime / setup.numWritesCompression));
                double nextWriteTime =
                        (currNumWrites + 1) * setup.compressionTime / setup.numWritesCompression + setup.dryingTime;

                if (nextTime + timeControl.GetMinTimeStep() > nextWriteTime)
                {
                    timestepHelper.prevTimestep = timestep;
                    timestepHelper.isTimestepAdjusted = true;
                    timestep = nextWriteTime - currTime +
                               1e-10; // +1e-10 to avoid floating point problems (minimal timestep)
                }
            }

            return timestep;
        }
        else
        {
            timeControl.RestorePreviousTime();
            return timeControl.GetTimeStep() * 0.5;
        }
    });

    NM.GetTimeControl().AdjustTimestep(0, 1, true);


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

    // Solve drrying cyclus
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.SetNumProcessors(4);
    S.CalculateMaximumIndependentSets();

    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, tolDisp);
    NM.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, tolNLES);
    NM.SetToleranceResidual(Node::eDof::RELATIVEHUMIDITY, tolRH);
    NM.SetToleranceResidual(Node::eDof::WATERVOLUMEFRACTION, tolWF);

    NM.AddCalculationStep({Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION});
    NM.AddCalculationStep({Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN});
    NM.SetMaxTimeStep(setup.dryingTime / setup.numWritesDrying);
    NM.SetMinTimeStep(1e-4);
    NM.SetTimeStep(setup.initialTimestepDrying);
    NM.SetPerformLineSearch(false);
    NM.SetMaxNumIterations(setup.dryingMaxIter);
    NM.PostProcessing().SetMinTimeStepPlot(setup.dryingTime / setup.numWritesDrying);
    NM.PostProcessing().SetResultDirectory(simulationName, true);
    NM.Solve(setup.dryingTime);


    NM.ClearActiveDofCalculationSteps();
    NM.AddCalculationStep({Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN});
    NM.SetMaxTimeStep(setup.maxCompressionDeltaTime);
    NM.SetMinTimeStep(setup.maxCompressionDeltaTime / 10000.);
    NM.SetTimeStep(setup.maxCompressionDeltaTime);
    NM.PostProcessing().SetMinTimeStepPlot(setup.compressionTime / setup.numWritesCompression);
    NM.Solve(setup.dryingTime + setup.compressionTime);
    outputFile.close();
}


// Main
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int NumArgs, char* Args[])
{

    try
    {
        if (NumArgs > 1)
        {
            std::cout << "Number of arguments: " << NumArgs - 1 << std::endl;
            for (int i = 1; i < NumArgs; ++i)
                std::cout << "Argument " << i << ": " << Args[i] << std::endl;

            SimulationSetup setup;
            setup.RH_env = ArgToDouble(Args[1]);

            if (NumArgs == 2)
            {
                std::cout << std::endl << "Starting simulation with the following parameters:" << std::endl;
                std::cout << "Environmental relative humidity = " << setup.RH_env << std::endl;
            }
        }
        else
        {
            SimulationSetup setup;
            setup.RH_env = 0.4;
            setup.dryingTime = 8.;
            setup.dryingTimeRamp = 0.5;
            setup.initialTimestepDrying = 0.1;
            setup.numWritesDrying = 20;

            //            setup.RH_env = 1.0;
            //            setup.dryingTime = 10.;
            //            setup.dryingTimeRamp = 10.;
            //            setup.numWritesDrying = 1;
            DS(setup, true, "DamageShrinkage3d");
        }
        return 0;
    }
    catch (Exception e)
    {

        std::cout << "NUTO exception occured:" << std::endl << e.what() << std::endl;
    }
    catch (std::exception e)
    {

        std::cout << "STD exception occured:" << std::endl << e.what() << std::endl;
    }
    catch (...)
    {

        std::cout << "UNKNOWN exception occured:" << std::endl;
    }
}
