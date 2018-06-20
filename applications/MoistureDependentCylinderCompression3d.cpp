#include "nuto/math/EigenSparseSolve.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"


#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportBoundary.h"
#include "integrands/MoistureTransportCoefficients.h"
#include "structures/MultiPhysicsStructureNew.h"

#include <fstream>
#include <iostream>

using namespace NuTo;
using namespace NuTo::Constraint;
using namespace NuTo::Integrands;
using namespace NuTo::Laws;
using namespace NuTo::Visualize;

class OutputFile
{
public:
    OutputFile(std::string fileName)
    {
        mFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
        if (!mFile.is_open())
            throw Exception(__PRETTY_FUNCTION__, "Can't open output file \"" + fileName + "\" !");
    }

    ~OutputFile()
    {
        mFile.close();
    }

    template <typename T>
    OutputFile& Write(const T& data)
    {
        mFile << data;
        mFile.flush();
        return *this;
    }

    template <typename T, typename... TArgs>
    OutputFile& Write(const T& data, const TArgs&... args)
    {
        mFile << data;
        Write(args...);
        return *this;
    }

private:
    std::ofstream mFile;
};

double TimeDependentDisplacement(double t)
{
    return -0.1 * t;
}

void SaveStressStrain(Group<CellInterface>& groupVolumeCellsTotal, GradientDamage<3>& integrandVolume,
                      DofType dofDisplacements, double t)
{
    static OutputFile file("MDCCStressStrain.dat");
    double V{0};
    Eigen::VectorXd Stress{Eigen::VectorXd::Zero(6)};
    for (CellInterface& cell : groupVolumeCellsTotal)
    {
        V += cell.Integrate([&](const CellIpData& cipd) { return 1.0; });
        auto test = cell.Integrate([&](const CellIpData& cipd) {
            DofVector<double> d;
            d[dofDisplacements] = integrandVolume.mLinearElasticDamage.Stress(
                    cipd.Apply(dofDisplacements, Nabla::Strain()),
                    integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)));
            return d;
        });
        Stress += test[dofDisplacements];
    }
    double AxialStress{Stress[2] / V};
    double AxialStrain{TimeDependentDisplacement(t) / 300};

    file.Write(AxialStress, ", ", AxialStrain, "\n");

    std::cout << "Stress: " << AxialStress << std::endl;
    std::cout << "Strain: " << AxialStrain << std::endl;
}


int main(int argc, char* argv[])
{
    std::cout << "Load mesh..." << std::endl;
    auto meshGmsh = MeshGmsh{"Cylinder.msh"};
    auto& mesh = meshGmsh.GetMeshFEM();
    MultiPhysicsStructure MPS{mesh};

    // Get element groups
    const auto& groupVolumeElementsTotal = meshGmsh.GetPhysicalGroup("Volume");
    const auto& groupSurfaceElementsUpperFace = meshGmsh.GetPhysicalGroup("UpperFace");
    const auto& groupSurfaceElementsLowerFace = meshGmsh.GetPhysicalGroup("LowerFace");
    const auto& groupSurfaceElementsSideFace = meshGmsh.GetPhysicalGroup("SideFace");
    const auto groupSurfaceElementsTotal =
            Unite(groupSurfaceElementsLowerFace, groupSurfaceElementsSideFace, groupSurfaceElementsUpperFace);

    // Create dof types
    const auto& dofDisplacements = MPS.AddDofType("displacements", 3);
    ScalarDofType dofNonLocal{"non local equivalent strain"};
    const auto& dofWaterVolumeFraction = MPS.AddDofType("water volume fraction", 1);
    const auto& dofRelativeHumidity = MPS.AddDofType("relative humidity", 1);
    std::vector<DofType> dofs_MT{dofWaterVolumeFraction, dofRelativeHumidity};

    // Create interpolation
    auto& interpolationTetrahedronQuadratic = mesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    auto& interpolationTriangleQuadratic = mesh.CreateInterpolation(InterpolationTriangleQuadratic());

    // Create integrations
    const auto& integrationTetrahedron3 = MPS.AddIntegrationType(IntegrationTypeTetrahedron(5));
    const auto& integrationTriangle5 = MPS.AddIntegrationType(IntegrationTypeTriangle(5));


    // Create dofs --------------------------------------------------------------------------------

    std::cout << "Create dofs..." << std::endl;

    AddDofInterpolation(&mesh, dofDisplacements, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofNonLocal, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupVolumeElementsTotal, interpolationTetrahedronQuadratic);
    AddDofInterpolation(&mesh, dofWaterVolumeFraction, groupSurfaceElementsTotal, interpolationTriangleQuadratic);
    AddDofInterpolation(&mesh, dofRelativeHumidity, groupSurfaceElementsTotal, interpolationTriangleQuadratic);


    // Set time derivatives
    MPS.SetNumTimeDerivatives(2);


    // Set nodal values ---------------------------------------------------------------------------

    std::cout << "Set nodal values..." << std::endl;

    //    MPS.SetNodalValues(dofRelativeHumidity, {1.0});
    //    MPS.SetNodalValues(dofWaterVolumeFraction, {0.1});
    std::cout << "Nodal values set..." << std::endl;

    // Create cells -------------------------------------------------------------------------------

    std::cout << "Create cells..." << std::endl;

    auto groupVolumeCellsTotal = MPS.CreateCells(groupVolumeElementsTotal, integrationTetrahedron3);
    auto groupSurfaceCellsTotal = MPS.CreateCells(groupSurfaceElementsTotal, integrationTriangle5);

    // Create integrands --------------------------------------------------------------------------

    std::cout << "Create integrants..." << std::endl;

    auto& integrandMoistureTransport = MPS.AddIntegrand(IntegrandWrapper<MoistureTransport, double>{
            MoistureTransport{dofWaterVolumeFraction, dofRelativeHumidity, MTCoefficientConstant{1.},
                              MTCoefficientConstant{0.01}, MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, 1.,
                              0.1, 0.2},
            &MoistureTransport::Gradient, &MoistureTransport::Stiffness, &MoistureTransport::Damping});

    auto& integrandMoistureTransportBoundary = MPS.AddIntegrand(IntegrandWrapper<MoistureTransportBoundary, double>{
            MoistureTransportBoundary{integrandMoistureTransport.GetIntegrand(), dofWaterVolumeFraction,
                                      dofRelativeHumidity, 2., 2., 0.4},
            &MoistureTransportBoundary::Gradient, &MoistureTransportBoundary::Stiffness});


    //    auto lawVolume = LinearElastic<3>{30.e9, 0.2};
    //    MomentumBalance<3> integrandVolume{dofDisplacements, lawVolume};

    auto material = Material::DefaultConcrete();
    material.c = 0.25;
    material.gf *= 0.25;
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
    constraints.Add(dofDisplacements, Direction(groupNodesTop, Eigen::Vector3d{0., 0., 1.}, TimeDependentDisplacement));


    MPS.RenumberDofs();

    // Setup time integration ---------------------------------------------------------------------

    std::cout << "Setup time integration scheme..." << std::endl;

    double t = 0.;
    double t_final = 40.;
    double dt_ini = t_final / 160.;
    double dt = dt_ini;

    // Newmark

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    // Quasi static

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
    QSS.mTolerance = 1e-4;


    // Visualization ------------------------------------------------------------------------------

    std::cout << "Setup visualization..." << std::endl;

    PostProcess pp("MDCC3dResults");
    pp.DefineVisualizer("Volume", groupVolumeCellsTotal, AverageHandler());
    pp.Add("Volume", dofDisplacements);
    pp.Add("Volume", dofNonLocal);
    //    pp.Add("Volume", dofWaterVolumeFraction);
    //    pp.Add("Volume", dofRelativeHumidity);
    pp.Add("Volume",
           [&](const CellIpData& cipd) { return integrandVolume.mDamageLaw.Damage(integrandVolume.Kappa(cipd)); },
           "Damage");

    std::cout << "Visualize first timestep..." << std::endl;
    pp.Plot(t, false);
    SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t);


    // Solve --------------------------------------------------------------------------------------

    std::cout << "Extract nodal values..." << std::endl;

    Eigen::VectorXd d_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 0);
    Eigen::VectorXd v_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 1);
    Eigen::VectorXd a_MT = MPS.ExtractDofs({dofWaterVolumeFraction, dofRelativeHumidity}, 2);
    Eigen::VectorXd delta_d_MT = Eigen::VectorXd::Zero(d_MT.rows());
    EigenSparseSolver solver("MumpsLU");


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


        //        Eigen::VectorXd d_it = d_MT + delta_d_MT;
        //        Eigen::VectorXd v_it =
        //                gamma / (dt * beta) * delta_d_MT + (1 - gamma / beta) * v_MT + dt * (1 - gamma / (2 * beta)) *
        //                a_MT;
        //        Eigen::VectorXd a_it = 1 / (dt * dt * beta) * delta_d_MT - 1 / (dt * beta) * v_MT - (1 / (2 * beta) -
        //        1) * a_MT;

        //        MPS.MergeDofs(d_it, dofs_MT, 0);
        //        MPS.MergeDofs(v_it, dofs_MT, 1);
        //        MPS.MergeDofs(a_it, dofs_MT, 2);


        //        Eigen::VectorXd grad = MPS.AssembleResidual(dofs_MT,
        //                                                    {{&groupVolumeCellsTotal, &integrandMoistureTransport},
        //                                                     {&groupSurfaceCellsTotal,
        //                                                     &integrandMoistureTransportBoundary}});

        //        int iteration = 0;
        //        while (grad.lpNorm<Eigen::Infinity>() > 10e-13)
        //        {

        //            std::cout << "Assemble moisture problem..." << std::endl;
        //            ++iteration;
        //            if (iteration > 20)
        //                throw Exception(__PRETTY_FUNCTION__, "No convergence");

        //            Eigen::SparseMatrix<double> K =
        //                    MPS.AssembleStiffness(dofs_MT,
        //                                          {{&groupVolumeCellsTotal, &integrandMoistureTransport},
        //                                           {&groupSurfaceCellsTotal, &integrandMoistureTransportBoundary}});

        //            Eigen::SparseMatrix<double> D =
        //                    MPS.AssembleDamping(dofs_MT, {{&groupVolumeCellsTotal, &integrandMoistureTransport}});
        //            Eigen::SparseMatrix<double> H = gamma / (dt * beta) * D + K;

        //            std::cout << "Solve moisture problem..." << std::endl;
        //            delta_d_MT = solver.Solve(H, -grad);

        //            d_it += delta_d_MT;
        //            v_it += gamma / (dt * beta) * delta_d_MT;
        //            a_it += 1 / (dt * dt * beta) * delta_d_MT;

        //            MPS.MergeDofs(d_it, dofs_MT, 0);
        //            MPS.MergeDofs(v_it, dofs_MT, 1);
        //            MPS.MergeDofs(a_it, dofs_MT, 2);

        //            grad = MPS.AssembleResidual(dofs_MT,
        //                                        {{&groupVolumeCellsTotal, &integrandMoistureTransport},
        //                                         {&groupSurfaceCellsTotal, &integrandMoistureTransportBoundary}});

        //            std::cout << "Time: " << t << std::endl
        //                      << "Iteration: " << iteration << std::endl
        //                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
        //                      << std::endl;
        //        }
        //        d_MT = d_it;
        //        v_MT = v_it;
        //        a_MT = a_it;


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
        if (dt < dt_ini && num_converges > 5)
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

        SaveStressStrain(groupVolumeCellsTotal, integrandVolume, dofDisplacements, t);
    }
}
