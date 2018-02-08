#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/elements/ElementBase.h"
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
    double dryingTimeRamp = 0.5;
    double initialTimestepDrying = 0.1;
    int dryingMaxIter = 10;

    // Compression Test
    double maxCompression = 2;
    double compressionTime = 1.;
    double maxCompressionDeltaTime = 0.0025;

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


int ArgToInt(std::string Arg)
{
    try
    {
        int number = std::stoi(Arg);
        return number;
    }
    catch (std::invalid_argument e)
    {
        throw Exception(__PRETTY_FUNCTION__, "Could not cast \"" + Arg + "\" into an integer!");
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
    if (1 < setup.RH_env)
        throw Exception(__PRETTY_FUNCTION__, "Relative humidity must be in the range [0,1]");

    if (0. >= setup.maxCompressionDeltaTime)
        throw Exception(__PRETTY_FUNCTION__, "Delta t of compression period must be > 0!");
    setup.maxCompression = std::abs(setup.maxCompression);


    // Adjust test name
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    simulationName = simulationName + "_RH" + std::to_string(setup.RH_env) + "_T" + std::to_string(setup.dryingTime);

    // Structure and TimeIntegration
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Structure S(3);
    NewmarkDirect NM(&S);

    S.SetShowTime(false);
    S.SetNumTimeDerivatives(1);
    S.SetVerboseLevel(0);

    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto meshData = S.ImportFromGmsh("ARAMIS3dgrob.msh");

    int elementGroupMatrix_Id = meshData[0].first;
    int interpolationTypeMatrix_Id = meshData[0].second;
    int elementGroupCylinder_Id = meshData[1].first;
    int interpolationTypeCylinder_Id = meshData[1].second;
    // auto prism = NuTo::MeshCompanion::ElementPrismsCreate(S, elementGroupMatrix_Id, elementGroupCylinder_Id, 0.1);

    std::vector<int> allNodes_Ids;
    S.NodeGroupGetMembers(S.GroupGetNodesTotal(), allNodes_Ids);


    // Constitutive laws
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // mechanics
    int lawLE_Id = S.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    S.ConstitutiveLawSetParameterDouble(lawLE_Id, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus * 2);
    S.ConstitutiveLawSetParameterDouble(lawLE_Id, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);

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
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH, 5.0e3);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV, 5.0e9);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DENSITY_WATER, 999.97);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH, 2.5e-1);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV, 1.0e5);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_RH, 1.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::DIFFUSION_EXPONENT_WV, 5.0);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION,
                                        0.56);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION,
                                        0.26);
    S.ConstitutiveLawSetParameterDouble(lawMT_Id, eConstitutiveParameter::MASS_EXCHANGE_RATE, 2.e-3);
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

    S.ElementGroupSetConstitutiveLaw(elementGroupMatrix_Id, lawAO_Id);
    S.ElementGroupSetConstitutiveLaw(elementGroupCylinder_Id, lawLE_Id);


    // Integration type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S.InterpolationTypeSetIntegrationType(interpolationTypeMatrix_Id,
                                          eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);
    S.InterpolationTypeSetIntegrationType(interpolationTypeCylinder_Id,
                                          eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);

    // Interpolation type
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // matrix
    S.InterpolationTypeAdd(interpolationTypeMatrix_Id, Node::eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationTypeMatrix_Id, Node::eDof::NONLOCALEQSTRAIN,
                           Interpolation::eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationTypeMatrix_Id, Node::eDof::RELATIVEHUMIDITY,
                           Interpolation::eTypeOrder::EQUIDISTANT2);
    S.InterpolationTypeAdd(interpolationTypeMatrix_Id, Node::eDof::WATERVOLUMEFRACTION,
                           Interpolation::eTypeOrder::EQUIDISTANT2);


    S.InterpolationTypeAdd(interpolationTypeCylinder_Id, Node::eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);

    // MeshCompanion::ElementPrismsCreate()
    S.ElementTotalConvertToInterpolationType();

    // Create necessary groups before mesh mapping
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto& bottomNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::Y, 0. - 1e-6, 0. + 1e-6);
    auto& topNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::Y, 200 - 1e-6, 200 + 1e+6);
    int topElementsGroup_Id = S.GroupCreateElementGroup();
    int bottomElementsGroup_Id = S.GroupCreateElementGroup();
    S.GroupAddElementsFromNodes(topElementsGroup_Id, S.GroupGetId(&topNodesGroup), false);
    S.GroupAddElementsFromNodes(bottomElementsGroup_Id, S.GroupGetId(&bottomNodesGroup), false);


    // Boundary Elements
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int BoundaryNodesGroup = S.GroupCreateNodeGroup();
    S.GroupAddNodeFunction(BoundaryNodesGroup, [](NodeBase* nodePtr) -> bool {
        auto coords = nodePtr->Get(Node::eDof::COORDINATES);
        if (nodePtr->IsDof(NuTo::Node::eDof::RELATIVEHUMIDITY) || nodePtr->IsDof(NuTo::Node::eDof::WATERVOLUMEFRACTION))
            if (coords[0] >= 200 - 1e-3 || coords[0] <= 0 + 1e-3 || coords[1] >= 200 - 1e-3 || coords[1] <= 0 + 1e-3 ||
                coords[2] >= 100 - 1e-3 || coords[2] <= 0 + 1e-3)
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
    std::vector<int> allElementIDsVector;

    S.ElementGroupGetMembers(allElementsGroup_Id, allElementIDsVector);

    std::vector<int> cylinderElementIDsVector;
    S.ElementGroupGetMembers(elementGroupCylinder_Id, cylinderElementIDsVector);

    for (unsigned int i = 0; i < allElementIDsVector.size(); i++)
    {

        for (int theIP = 0; theIP < S.ElementGetElementPtr(allElementIDsVector[i])->GetNumIntegrationPoints(); theIP++)
        {

            // TODO: Skip elements of cylinders!
            if (S.ElementGetElementPtr(i)->GetConstitutiveLaw(theIP).GetType() ==
                eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS)
                break;
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

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Y}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::Y},
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
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(bottomNodesGroup, {eDirection::Z}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::X}));
    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(topNodesGroup, {eDirection::Z}));

    // Postprocessing for compression test
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    double topElementsVolume = S.ElementGroupGetVolume(topElementsGroup_Id);
    double bottomElementsVolume = S.ElementGroupGetVolume(bottomElementsGroup_Id);
    //    S.ElementGroupGetAverageStress(topElementsGroup_Id, S.ElementGroupGetVolume(topElementsGroup_Id),
    //    StressVector);


    std::ofstream outputFile;
    outputFile.open(simulationName + ".txt");
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

    int visualizeGroup = elementGroupMatrix_Id;
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
    NM.SetMinTimeStep(setup.maxCompressionDeltaTime / 1000.);
    NM.SetTimeStep(setup.maxCompressionDeltaTime);
    NM.PostProcessing().SetMinTimeStepPlot(setup.compressionTime / setup.numWritesCompression);
    NM.Solve(setup.dryingTime + setup.compressionTime);
    outputFile.close();
}


// Main
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int NumArgs, char* Args[])
{
    constexpr int neededNumArgs = 6;
    try
    {
        if (NumArgs > 1)
        {
            SimulationSetup setup;
            if (Args[1][0] == '-')
            {
                std::string command(Args[1]);
                if (command.compare("-h") == 0 || command.compare("-help") == 0)
                {
                    std::cout << "Necessary number of arguments:" << neededNumArgs << std::endl << std::endl;
                    std::cout << "Argument 1: Environmental relative humidity [0,1]" << std::endl;
                    std::cout << "Argument 2: drying Time (> " << setup.dryingTimeRamp << ")" << std::endl;
                    std::cout << "Argument 3: Number of postprocessings during drying period" << std::endl;
                    std::cout << "Argument 4: Minimum number of timesteps during compression period (should be > 20)"
                              << std::endl;
                    std::cout << "Argument 5: Number of postprocessings during compression period" << std::endl;
                    return 0;
                }
                else
                {
                    std::cout << "Unknown command: \"" << Args[1] << "\"" << std::endl;
                    return 0;
                }
            }

            std::cout << "Number of arguments: " << NumArgs - 1 << std::endl;
            for (int i = 1; i < NumArgs; ++i)
                std::cout << "Argument " << i << ": " << Args[i] << std::endl;


            if (NumArgs == neededNumArgs)
            {
                setup.RH_env = ArgToDouble(Args[1]);
                setup.dryingTime = ArgToDouble(Args[2]);
                setup.numWritesDrying = ArgToInt(Args[3]);
                setup.maxCompressionDeltaTime = setup.maxCompression / std::round(ArgToDouble(Args[4]));
                setup.numWritesCompression = ArgToInt(Args[5]);

                std::cout << std::endl << "Starting simulation with the following parameters:" << std::endl;
                std::cout << "Environmental relative humidity = " << setup.RH_env << " [0,1]" << std::endl;
                std::cout << "Length of the drying period = " << setup.dryingTime << " (> " << setup.dryingTimeRamp
                          << ")" << std::endl;
                std::cout << "Number of postprocessings during drying period = " << setup.numWritesDrying << std::endl;
                std::cout << "Maximum timestep during compression period = " << setup.maxCompressionDeltaTime
                          << " (should be < " << setup.compressionTime / 20. << ")" << std::endl;
                std::cout << "Number of postprocessings during compression period = " << setup.numWritesCompression
                          << std::endl;
                DS(setup, true, "DS");
            }
            else
            {
                std::cout << std::endl
                          << "Incorrect number of parameters: " << NumArgs - 1 << std::endl
                          << "Necessary number of arguments:" << neededNumArgs - 1 << std::endl
                          << "Use -h to get parameter Info ... Closing program!" << std::endl;
            }
        }
        else
        {
            SimulationSetup setup;
            setup.RH_env = 0.4;
            setup.dryingTime = 50.;
            setup.dryingTimeRamp = 0.5;
            setup.initialTimestepDrying = 0.1;
            setup.numWritesDrying = 25;

            //            setup.RH_env = 1.0;
            //            setup.dryingTime = 10.;
            //            setup.dryingTimeRamp = 10.;
            //            // setup.initialTimestepDrying = 5.;
            //            setup.numWritesDrying = 10;
            //            setup.maxCompressionDeltaTime = 0.01;

            DS(setup, true, "Aramis3d");
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
