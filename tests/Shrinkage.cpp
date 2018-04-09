#include <boost/test/unit_test.hpp>

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/dofs/DofNumbering.h"
#include "nuto/mechanics/dofs/GlobalDofVector.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/tools/NodalValueMerger.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PostProcess.h"


#include "integrands/MoistureTransportCoefficients.h"
#include "integrands/MoistureTransport.h"
#include "integrands/Shrinkage.h"

#include <boost/ptr_container/ptr_vector.hpp>

using namespace NuTo;
using namespace NuTo::Integrands;
using namespace std::placeholders;
using MoistureTransportIntegrand = MoistureTransport<2, MTCConst<1>, MTCConst<1, 100>, MTCConst<0>, MTCConst<1, 10>>;

class ShrinkageTest
{
    DofType mDofDisp = {"Displacements", 3};
    DofType mDofRH = {"relativeHumidity", 1};
    DofType mDofWV = {"waterVolumeFraction", 1};
    DofInfo mDofInfo;
    Constraint::Constraints mConstraints;
    IntegrationTypeTriangle integrationType = {2};
    MeshFem mMesh;
    Group<CellInterface> mGrpCellsMT;
    boost::ptr_vector<CellInterface> mCellContainer;
    MoistureTransportIntegrand mMoistureTransport;

public:
    ShrinkageTest(double rho_w = 1., double rho_g_sat = 0.1, double PV = 0.2)
        : mMoistureTransport(mDofWV, mDofRH, rho_w, rho_g_sat, PV)
    {
    }

    std::vector<DofType> GetDofs()
    {
        return {mDofDisp, mDofWV, mDofRH};
    }

    const Group<CellInterface> GetMoistureTransportCells() const
    {
        return mGrpCellsMT;
    }

    MeshFem& GetMesh()
    {
        return mMesh;
    }


    GlobalDofVector GradientMoistureTransport()
    {
        CheckDofNumbering();
        return SimpleAssembler(mDofInfo).BuildVector(
                mGrpCellsMT, {mDofRH, mDofWV},
                std::bind(&MoistureTransportIntegrand::Gradient, mMoistureTransport, _1, 0.));
    }

    GlobalDofMatrixSparse Stiffness();
    GlobalDofMatrixSparse Damping();

    void CheckDofNumbering();
    GlobalDofVector CreateGlobalDofVector();
    void GetDofVector(GlobalDofVector& dofs, int instance = 0);

    void MergeDofs(Eigen::VectorXd d, Eigen::VectorXd v, Eigen::VectorXd a);
    void MergeDofVector(GlobalDofVector& dofs, int instance = 0);

    void CreateUnitMesh(int numX, int numY, const InterpolationSimple& interpolationDispArg,
                        const InterpolationSimple& interpolationWVArg, const InterpolationSimple& interpolationRHArg);

    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ExtractDofs();
};


GlobalDofMatrixSparse ShrinkageTest::Stiffness()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildMatrix(
            mGrpCellsMT, {mDofRH, mDofWV},
            std::bind(&MoistureTransportIntegrand::Stiffness, mMoistureTransport, _1, 0.));
}

GlobalDofMatrixSparse ShrinkageTest::Damping()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildMatrix(
            mGrpCellsMT, {mDofRH, mDofWV}, std::bind(&MoistureTransportIntegrand::Damping, mMoistureTransport, _1, 0.));
}

void ShrinkageTest::CheckDofNumbering()
{
    if (!(mDofInfo.numDependentDofs.Has(mDofDisp) && mDofInfo.numDependentDofs.Has(mDofWV) &&
          mDofInfo.numDependentDofs.Has(mDofRH)))
        throw Exception(__PRETTY_FUNCTION__, "Please do a dof numbering before trying to assemble something!");
}

GlobalDofVector ShrinkageTest::CreateGlobalDofVector()
{
    GlobalDofVector d;
    //    int nDofDisp = (mMesh.NodesTotal(mDofDisp)).Size() * mDofDisp.GetNum();
    //    int nDofDisp_dep = mConstraints.GetNumEquations(mDofDisp);
    int nDofWV = (mMesh.NodesTotal(mDofWV)).Size() * mDofWV.GetNum();
    int nDofWV_dep = mConstraints.GetNumEquations(mDofWV);
    int nDofRH = (mMesh.NodesTotal(mDofRH)).Size() * mDofRH.GetNum();
    int nDofRH_dep = mConstraints.GetNumEquations(mDofRH);

    //    d.J[mDofDisp] = Eigen::VectorXd(nDofDisp - nDofDisp_dep);
    //    d.K[mDofDisp] = Eigen::VectorXd(nDofDisp_dep);
    d.J[mDofWV] = Eigen::VectorXd(nDofWV - nDofWV_dep);
    d.K[mDofWV] = Eigen::VectorXd(nDofWV_dep);
    d.J[mDofRH] = Eigen::VectorXd(nDofRH - nDofRH_dep);
    d.K[mDofRH] = Eigen::VectorXd(nDofRH_dep);

    return d;
}

void ShrinkageTest::GetDofVector(GlobalDofVector& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Extract(&dofs, {mDofWV, mDofRH}, instance);
}

void ShrinkageTest::MergeDofs(Eigen::VectorXd d, Eigen::VectorXd v, Eigen::VectorXd a)
{
    GlobalDofVector dg{CreateGlobalDofVector()}, vg{CreateGlobalDofVector()}, ag{CreateGlobalDofVector()};
    FromEigen(d, {mDofWV, mDofRH}, &dg.J);
    FromEigen(v, {mDofWV, mDofRH}, &vg.J);
    FromEigen(a, {mDofWV, mDofRH}, &ag.J);
    MergeDofVector(dg);
    MergeDofVector(vg, 1);
    MergeDofVector(ag, 2);
}

void ShrinkageTest::MergeDofVector(GlobalDofVector& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Merge(dofs, dofs.J.DofTypes(), instance);
}

void ShrinkageTest::CreateUnitMesh(int numX, int numY, const InterpolationSimple& interpolationDispArg,
                                   const InterpolationSimple& interpolationWVArg,
                                   const InterpolationSimple& interpolationRHArg)
{
    mMesh = UnitMeshFem::CreateTriangles(numX, numY);

    // Interpolation
    const auto& interpolationDisp = mMesh.CreateInterpolation(interpolationDispArg);
    const auto& interpolationWV = mMesh.CreateInterpolation(interpolationWVArg);
    const auto& interpolationRH = mMesh.CreateInterpolation(interpolationRHArg);

    // Add Dofs
    AddDofInterpolation(&mMesh, mDofDisp, interpolationDisp);
    AddDofInterpolation(&mMesh, mDofWV, interpolationWV);
    AddDofInterpolation(&mMesh, mDofRH, interpolationRH);
    mMesh.AllocateDofInstances(mDofDisp, 3);
    mMesh.AllocateDofInstances(mDofWV, 3);
    mMesh.AllocateDofInstances(mDofRH, 3);


    //         dof numbering
    mDofInfo = DofNumbering::Build(mMesh.NodesTotal(mDofDisp), mDofDisp, mConstraints);
    mDofInfo.Merge(mDofWV, DofNumbering::Build(mMesh.NodesTotal(mDofWV), mDofWV, mConstraints));
    mDofInfo.Merge(mDofRH, DofNumbering::Build(mMesh.NodesTotal(mDofRH), mDofRH, mConstraints));

    // set nodal values to 100% RH and equilibium state
    for (NodeSimple& node : mMesh.NodesTotal(mDofWV))
    {
        node.SetValue(0, 0.10);
        node.SetValue(0, 0., 1);
        node.SetValue(0, 0., 2);
    }
    for (NodeSimple& node : mMesh.NodesTotal(mDofRH))
    {
        node.SetValue(0, 1.);
        node.SetValue(0, 0., 1);
        node.SetValue(0, 0., 2);
    }


    // create cells
    int cellId = 0;
    for (ElementCollection& element : mMesh.Elements)
    {
        mCellContainer.push_back(new Cell(element, integrationType, cellId++));
        mGrpCellsMT.Add(mCellContainer.back());
    }

    //    // create Boundary elements
    //    InterpolationSimple& boundaryInterpolation = mMesh.CreateInterpolation(InterpolationPoint());

    //    Eigen::VectorXd coordLeftBoundary = Eigen::VectorXd::Zero(1);
    //    Eigen::VectorXd coordRightBoundary = Eigen::VectorXd::Ones(1);

    //    ElementCollectionFem& leftElementCollection =
    //            mMesh.Elements.Add({{{mMesh.NodeAtCoordinate(coordLeftBoundary)}, boundaryInterpolation}});
    //    leftElementCollection.AddDofElement(
    //            dofWV, ElementFem({&mMesh.NodeAtCoordinate(coordLeftBoundary, dofWV)}, boundaryInterpolation));
    //    leftElementCollection.AddDofElement(
    //            dofRH, ElementFem({&mMesh.NodeAtCoordinate(coordLeftBoundary, dofRH)}, boundaryInterpolation));


    //    ElementCollectionFem& rightElementCollection =
    //            mMesh.Elements.Add({{{mMesh.NodeAtCoordinate(coordRightBoundary)}, boundaryInterpolation}});
    //    rightElementCollection.AddDofElement(
    //            dofWV, ElementFem({&mMesh.NodeAtCoordinate(coordRightBoundary, dofWV)}, boundaryInterpolation));
    //    rightElementCollection.AddDofElement(
    //            dofRH, ElementFem({&mMesh.NodeAtCoordinate(coordRightBoundary, dofRH)}, boundaryInterpolation));

    //    cellContainer.push_back(new CellPoint(leftElementCollection, cellId++));
    //    grpCellsMTBoundary.Add(cellContainer.back());
    //    cellContainer.push_back(new CellPoint(rightElementCollection, cellId++));
    //    grpCellsMTBoundary.Add(cellContainer.back());
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ShrinkageTest::ExtractDofs()
{
    GlobalDofVector dg{CreateGlobalDofVector()}, vg{CreateGlobalDofVector()}, ag{CreateGlobalDofVector()};
    GetDofVector(dg);
    GetDofVector(vg, 1);
    GetDofVector(ag, 2);
    return std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>(
            ToEigen(dg.J, {mDofWV, mDofRH}), ToEigen(vg.J, {mDofWV, mDofRH}), ToEigen(ag.J, {mDofWV, mDofRH}));
}

BOOST_AUTO_TEST_CASE(Construction)
{
    ShrinkageTest ST;
    ST.CreateUnitMesh(10, 10, InterpolationTriangleQuadratic(), InterpolationTriangleQuadratic(),
                      InterpolationTriangleQuadratic());

    Eigen::VectorXd d_MT, v_MT, a_MT;
    std::tie(d_MT, v_MT, a_MT) = ST.ExtractDofs();

    auto dofs = ST.GetDofs();
    const int numDofsDisp = ST.GetMesh().NodesTotal(dofs[0]).Size() * dofs[0].GetNum();
    const int numDofsWV = ST.GetMesh().NodesTotal(dofs[1]).Size() * dofs[1].GetNum();
    const int numDofsRH = ST.GetMesh().NodesTotal(dofs[2]).Size() * dofs[2].GetNum();


    for (int i = 0; i < numDofsWV; ++i)
    {
        d_MT[i] = 0.1;
    }
    for (int i = 0; i < numDofsRH; ++i)
    {
        d_MT[i + numDofsWV] = 0.4 + 0.6 * i / numDofsRH;
    }

    ST.MergeDofs(d_MT, v_MT, a_MT);

    Visualize::PostProcess pp("ShrinkageResults");
    pp.DefineVisualizer("Moisture Transport", ST.GetMoistureTransportCells(), Visualize::AverageHandler());
    pp.Add("Moisture Transport", dofs[1]);
    pp.Add("Moisture Transport", dofs[2]);
    pp.Plot(0., false);
    // visu.DefineCellData("Tensor");

    constexpr double gamma = 1. / 2.;
    constexpr double beta = 1. / 4.;


    EigenSparseSolver solver("EigenSparseLU");
    Eigen::VectorXd delta_d = Eigen::VectorXd::Zero(d_MT.rows());

    double t = 0.;
    double delta_t = 0.01;
    double t_final = 0.1;
    while (t < t_final)
    {
        t += delta_t;
        int iteration = 0;

        //        //        Eigen::VectorXd gradBoundary = ToEigen(MTT.GradientBoundary().J, MTT.GetDofs());
        //        //        Eigen::MatrixXd stiffBoundary = ToEigen(MTT.StiffnessBoundary().JJ, MTT.GetDofs());
        //        //        std::cout << stiffBoundary << std::endl;

        Eigen::VectorXd d_it = d_MT + delta_d;
        Eigen::VectorXd v_it = gamma / (delta_t * beta) * delta_d + (1 - gamma / beta) * v_MT +
                               delta_t * (1 - gamma / (2 * beta)) * a_MT;
        Eigen::VectorXd a_it =
                1 / (delta_t * delta_t * beta) * delta_d - 1 / (delta_t * beta) * v_MT - (1 / (2 * beta) - 1) * a_MT;
        ST.MergeDofs(d_it, v_it, a_it);


        Eigen::VectorXd grad = ToEigen(ST.GradientMoistureTransport().J, {dofs[1], dofs[2]});

        while (grad.lpNorm<Eigen::Infinity>() > 10e-12)
        {

            ++iteration;
            if (iteration > 20)
                throw Exception(__PRETTY_FUNCTION__, "No convergence");

            Eigen::SparseMatrix<double> K = ToEigen(ST.Stiffness().JJ, {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> D = ToEigen(ST.Damping().JJ, {dofs[1], dofs[2]});
            Eigen::SparseMatrix<double> H = gamma / (delta_t * beta) * D + K;

            delta_d = solver.Solve(H, -grad);

            d_it += delta_d;
            v_it += gamma / (delta_t * beta) * delta_d;
            a_it += 1 / (delta_t * delta_t * beta) * delta_d;

            ST.MergeDofs(d_it, v_it, a_it);

            grad = ToEigen(ST.GradientMoistureTransport().J, {dofs[1], dofs[2]});
            std::cout << "Time: " << t << std::endl
                      << "Iteration: " << iteration << std::endl
                      << "Max residual: " << grad.lpNorm<Eigen::Infinity>() << std::endl
                      << std::endl;
        }
        d_MT = d_it;
        v_MT = v_it;
        a_MT = a_it;
        pp.Plot(t, false);
    }
}
