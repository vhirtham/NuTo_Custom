#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
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
    auto meshGmsh = MeshGmsh{"ARAMIS2d.msh"};
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};

    // Get element groups
    const auto& groupSurfaceElementsMatrix = meshGmsh.GetPhysicalGroup("matrix");
    const auto& groupSurfaceElementsGranite = meshGmsh.GetPhysicalGroup("granite");
    const auto& groupSurfaceElementsTotal = Unite(groupSurfaceElementsMatrix, groupSurfaceElementsGranite);

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("displacements", 2);
    ScalarDofType dofNonLocal{"non local equivalent strain"};

    // Create interpolation
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTriangle5 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofDisplacements, groupSurfaceElementsTotal, interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofNonLocal, groupSurfaceElementsMatrix, interpolationTriangleQuadratic);


    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupSurfaceCellsGranite = MPS.CreateCells(groupSurfaceElementsGranite, integrationTriangle5);
    auto groupSurfaceCellsMatrix = MPS.CreateCells(groupSurfaceElementsMatrix, integrationTriangle5);
    auto groupSurfaceCellsTotal = Unite(groupSurfaceCellsGranite, groupSurfaceCellsMatrix);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;

    //    auto lawMatrix = LinearElastic<2>{30.e9, 0.2};
    //    MomentumBalance<2> integrandMatrix{dofDisplacements, lawMatrix};

    auto material = Material::DefaultConcrete();
    material.c = 0.25;
    GradientDamage<2> integrandMatrix{dofDisplacements, dofNonLocal, material};
    int nIp = integrationTriangle5.GetNumIntegrationPoints();
    integrandMatrix.mKappas.resize(groupSurfaceCellsTotal.Size(), nIp);

    auto lawGranite = LinearElastic<2>{60000., 0.2};
    MomentumBalance<2> integrandGranite{dofDisplacements, lawGranite};


    // Define constraints -------------------------------------------------------------------------

    Constraints& constraints = MPS.GetConstraints();

    auto& nodeBottomLeft = mesh.NodeAtCoordinate(Eigen::Vector2d{0., 0.}, dofDisplacements);
    auto groupNodesBottom = mesh.NodesAtAxis(eDirection::Y, dofDisplacements);
    auto groupNodesTop = mesh.NodesAtAxis(eDirection::Y, dofDisplacements, 200);

    constraints.Add(dofDisplacements, Direction(nodeBottomLeft, Eigen::Vector2d{1., 0.}));
    constraints.Add(dofDisplacements, Direction(groupNodesBottom, Eigen::Vector2d{0., 1.}));
    constraints.Add(dofDisplacements,
                    Direction(groupNodesTop, Eigen::Vector2d{0., 1.}, [](double t) { return 0.01 * t; }));


    // Setup time integration ---------------------------------------------------------------------

    std::cout << "Setup time integration scheme..." << std::endl;

    double t = 0.;
    double t_final = 40.;
    double dt_ini = t_final / 40.;
    double dt = dt_ini;


    TimeDependentProblem TDP{&mesh};
    TDP.AddGradientFunction(groupSurfaceCellsMatrix, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMatrix.Gradient(cipd);
    });
    TDP.AddHessian0Function(groupSurfaceCellsMatrix, [&](const CellIpData& cipd, double t, double dt) {
        return integrandMatrix.Hessian0(cipd);
    });
    TDP.AddUpdateFunction(groupSurfaceCellsMatrix,
                          [&](const CellIpData& cipd, double t, double dt) { return integrandMatrix.Update(cipd); });

    TDP.AddGradientFunction(groupSurfaceCellsGranite, [&](const CellIpData& cipd, double t, double dt) {
        return integrandGranite.Gradient(cipd, dt);
    });
    TDP.AddHessian0Function(groupSurfaceCellsGranite, [&](const CellIpData& cipd, double t, double dt) {
        return integrandGranite.Hessian0(cipd, dt);
    });

    QuasistaticSolver QSS{TDP, {dofDisplacements, dofNonLocal}};
    QSS.SetGlobalTime(t);
    QSS.SetConstraints(constraints);
    QSS.mTolerance = 1e-6;


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    PostProcess pp("GradientDamageAramis2dResults");
    pp.DefineVisualizer("Matrix", groupSurfaceCellsMatrix, AverageHandler());
    pp.Add("Matrix", dofDisplacements);
    pp.Add("Matrix", dofNonLocal);
    pp.Add("Matrix",
           [&](const CellIpData& cipd) { return integrandMatrix.mDamageLaw.Damage(integrandMatrix.Kappa(cipd)); },
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
