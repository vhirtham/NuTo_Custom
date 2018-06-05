#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"


#include "structures/MultiPhysicsStructureNew.h"

#include <iostream>

using namespace NuTo;
using namespace NuTo::Constraint;
using namespace NuTo::Integrands;
using namespace NuTo::Laws;
using namespace NuTo::Visualize;

int main(int argc, char* argv[])
{
    std::cout << "Load mesh..." << std::endl;
    auto meshGmsh = MeshGmsh{"Cylinder.msh"};
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};

    // Get element groups
    const auto& groupVolumeElementsTotal = meshGmsh.GetPhysicalGroup("Volume");

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("displacements", 3);
    ScalarDofType dofNonLocal{"non local equivalent strain"};

    // Create interpolation
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());

    // Create integrations
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofNonLocal, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);


    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsTotal = MPS.CreateCells(groupVolumeElementsTotal, integrationTetrahedron3);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;

    //    auto lawVolume = LinearElastic<3>{30.e9, 0.2};
    //    MomentumBalance<3> integrandVolume{dofDisplacements, lawVolume};

    auto material = Material::DefaultConcrete();
    material.c = 0.25;
    GradientDamage<3> integrandVolume{dofDisplacements, dofNonLocal, material};
    int nIp = integrationTetrahedron3.GetNumIntegrationPoints();
    integrandVolume.mKappas.resize(groupVolumeCellsTotal.Size(), nIp);


    // Define constraints -------------------------------------------------------------------------

    Constraints& constraints = MPS.GetConstraints();

    //    auto& nodeBottomRight = mesh.NodeAtCoordinate(Eigen::Vector3d{50., 0., 0.}, dofDisplacements, 1e-2);
    //    auto& nodeBottomLeft = mesh.NodeAtCoordinate(Eigen::Vector3d{-50., 0., 0.}, dofDisplacements);
    //    auto& nodeBottomFront = mesh.NodeAtCoordinate(Eigen::Vector3d{0., 50., 0.}, dofDisplacements);
    //    auto& nodeBottomBack = mesh.NodeAtCoordinate(Eigen::Vector3d{0., -50., 0.}, dofDisplacements);
    auto groupNodesBottom = mesh.NodesAtAxis(eDirection::Z, dofDisplacements);
    auto groupNodesTop = mesh.NodesAtAxis(eDirection::Z, dofDisplacements, 300);

    //    constraints.Add(dofDisplacements,
    //                    {Direction(nodeBottomFront, Eigen::Vector3d{1., 0., 0.}),
    //                     Direction(nodeBottomRight, Eigen::Vector3d{0., 1., 0.}),
    //                     Direction(nodeBottomLeft, Eigen::Vector3d{0., 1., 0.}),
    //                     Direction(nodeBottomBack, Eigen::Vector3d{1., 0., 0.})});

    constraints.Add(dofDisplacements, Direction(groupNodesBottom, Eigen::Vector3d{1., 0., 0.}));
    constraints.Add(dofDisplacements, Direction(groupNodesBottom, Eigen::Vector3d{0., 1., 0.}));
    constraints.Add(dofDisplacements, Direction(groupNodesBottom, Eigen::Vector3d{0., 0., 1.}));

    constraints.Add(dofDisplacements, Direction(groupNodesTop, Eigen::Vector3d{1., 0., 0.}));
    constraints.Add(dofDisplacements, Direction(groupNodesTop, Eigen::Vector3d{0., 1., 0.}));
    constraints.Add(dofDisplacements,
                    Direction(groupNodesTop, Eigen::Vector3d{0., 0., 1.}, [](double t) { return -0.1 * t; }));


    // Setup time integration ---------------------------------------------------------------------

    std::cout << "Setup time integration scheme..." << std::endl;

    double t = 0.;
    double t_final = 40.;
    double dt_ini = t_final / 400.;
    double dt = dt_ini;


    TimeDependentProblem TDP{&mesh};
    TDP.AddGradientFunction(groupVolumeCellsTotal, [&](const CellIpData& cipd, double t, double dt) {
        return integrandVolume.Gradient(cipd);
    });
    TDP.AddHessian0Function(groupVolumeCellsTotal, [&](const CellIpData& cipd, double t, double dt) {
        return integrandVolume.Hessian0(cipd);
    });
    TDP.AddUpdateFunction(groupVolumeCellsTotal,
                          [&](const CellIpData& cipd, double t, double dt) { return integrandVolume.Update(cipd); });


    // QuasistaticSolver QSS{TDP, dofDisplacements};
    QuasistaticSolver QSS{TDP, {dofDisplacements, dofNonLocal}};
    QSS.SetGlobalTime(t);
    QSS.SetConstraints(constraints);
    QSS.mTolerance = 1e-5;


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    PostProcess pp("GradientDamageCylinder3dResults");
    pp.DefineVisualizer("Volume", groupVolumeCellsTotal, AverageHandler());
    pp.Add("Volume", dofDisplacements);
    pp.Add("Volume", dofNonLocal);
    pp.Add("Volume",
           [&](const CellIpData& cipd) { return integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)); },
           "Damage");

    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);


    // Solve --------------------------------------------------------------------------------------

    double next_plot = dt_ini;
    int num_converges = 0;

    std::cout << "Start time integration scheme..." << std::endl;
    while (t < t_final)
    {
        t += dt;
        std::cout << std::endl
                  << "Next timestep" << std::endl
                  << "-------------" << std::endl
                  << "t  : " << t << std::endl
                  << "dt : " << dt << std::endl;
        std::cout << "Assemble and solve mechanics problem..." << std::endl;
        try
        {
            QSS.DoStep(t, "MumpsLU");
        }
        catch (NewtonRaphson::NoConvergence e)
        {
            std::cout << "No convergence... reducing timestep" << std::endl;
            t -= dt;
            dt *= 0.5;
            num_converges = 0;
            if (dt < dt_ini * 1e-6)
                throw e;
            continue;
        }

        ++num_converges;
        if (dt < dt_ini && num_converges > 10)
        {
            std::cout << "Increasing timestep" << std::endl;
            dt *= 1.5;
            if (dt > dt_ini)
                dt = dt_ini;
            num_converges = 0;
        }

        if (t >= next_plot)
        {
            std::cout << std::endl << "Visualize..." << std::endl << std::endl;
            pp.Plot(t, false);
            next_plot += dt_ini;
        }
    }
}
