#include "mechanics/DirectionEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

#include <functional>

using namespace NuTo;

int main(int argc, char* argv[])
{
    Structure S(2);
    auto mshData = S.ImportFromGmsh("ARAMIS.msh");


    // create section
    double thickness = 20.0;
    auto section = SectionPlane::Create(thickness, true);
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
    auto GUpBound = S.GroupGetNodeCoordinateRange(eDirection::Y, 200.0, 200.0);
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
    S.AddVisualizationComponent(GIVis, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);


    // Solve

    NewmarkDirect NM(&S);
    NM.SetAutomaticTimeStepping(false);
    NM.SetTimeStep(1800.0);
    NM.SetPerformLineSearch(false);
    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-2);

    NM.PostProcessing().SetResultDirectory("ARAMISResults", true);
    NM.Solve(3600.0);


    return 0;
}
