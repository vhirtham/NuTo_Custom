#pragma once

// NuTo includes
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/dofs/DofNumbering.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/tools/NodalValueMerger.h"
#include "nuto/math/shapes/Pyramid.h"

// NuToCustom includes
#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportBoundary.h"

// Test includes
#include "CellPoint.h"
#include "InterpolationPoint.h"

// Other includes
#include <boost/ptr_container/ptr_vector.hpp>

using namespace NuTo;
using namespace NuTo::Integrands;
using namespace std::placeholders;


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
class MoistureTransportTest
{
public:
    MoistureTransportTest(double rho_w = 1.0, double rho_g_sat = 1.0, double PV = 0.25);

    void CreateUnitMesh(int numElements, const InterpolationSimple& interpolationWVArg,
                        const InterpolationSimple& interpolationRHArg);


    DofVector<double> Gradient();
    DofMatrixSparse<double> Stiffness();
    DofMatrixSparse<double> Damping();
    DofVector<double> GradientBoundary();
    DofMatrixSparse<double> StiffnessBoundary();


    DofVector<double> CreateDofVector();
    void GetDofVector(DofVector<double>& dofs, int instance = 0);
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ExtractDofs();
    void MergeDofVector(DofVector<double>& dofs, int instance = 0);
    void MergeDofs(Eigen::VectorXd d, Eigen::VectorXd v, Eigen::VectorXd a);

    std::vector<DofType> GetDofs()
    {
        return {dofWV, dofRH};
    }


    MeshFem& GetMesh() &
    {
        return mMesh;
    }


private:
    void CheckDofNumbering()
    {
        if (!(dofInfo.numDependentDofs.Has(dofWV) && dofInfo.numDependentDofs.Has(dofRH)))
            throw Exception(__PRETTY_FUNCTION__, "Please do a dof numbering before trying to assemble something!");
    }

    // members
    DofType dofRH = {"relativeHumidity", 1};
    DofType dofWV = {"waterVolumeFraction", 1};
    DofInfo dofInfo;
    Constraint::Constraints constraints;
    MeshFem mMesh;
    boost::ptr_vector<CellInterface> cellContainer;
    Group<CellInterface> grpCellsMT;
    Group<CellInterface> grpCellsMTBoundary;
    IntegrationTypeTensorProduct<1> integrationType = {3, eIntegrationMethod::GAUSS};
    MoistureTransport integrandMoistureTransport;
    MoistureTransportBoundary<TDim, TDCw, TDCg, TWVEq> mIntegrandMoistureTransportBoundary;
};


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MoistureTransportTest(double rho_w, double rho_g_sat, double PV)
    : integrandMoistureTransport(dofWV, dofRH, MTCoefficientConstant{1.}, MTCoefficientConstant{0.01},
                                 MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, rho_w, rho_g_sat, PV)
    , mIntegrandMoistureTransportBoundary(dofWV, dofRH)
{
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::CreateUnitMesh(int numElements,
                                                                          const InterpolationSimple& interpolationWVArg,
                                                                          const InterpolationSimple& interpolationRHArg)
{
    // mesh creation
    mMesh = UnitMeshFem::CreateLines(numElements);
    const auto& interpolationWV = mMesh.CreateInterpolation(interpolationWVArg);
    const auto& interpolationRH = mMesh.CreateInterpolation(interpolationRHArg);
    AddDofInterpolation(&mMesh, dofWV, interpolationWV);
    AddDofInterpolation(&mMesh, dofRH, interpolationRH);

    mMesh.AllocateDofInstances(dofWV, 3);
    mMesh.AllocateDofInstances(dofRH, 3);

    // dof numbering
    dofInfo = DofNumbering::Build(mMesh.NodesTotal(dofWV), dofWV, constraints);
    dofInfo.Merge(dofRH, DofNumbering::Build(mMesh.NodesTotal(dofRH), dofRH, constraints));

    // set nodal values to 100% RH and equilibium state
    for (NodeSimple& node : mMesh.NodesTotal(dofWV))
    {
        node.SetValue(0, TWVEq::value(0., 1., 0., 0.));
        node.SetValue(0, 0., 1);
        node.SetValue(0, 0., 2);
    }
    for (NodeSimple& node : mMesh.NodesTotal(dofRH))
    {
        node.SetValue(0, 1.);
        node.SetValue(0, 0., 1);
        node.SetValue(0, 0., 2);
    }


    // create cells
    int cellId = 0;
    for (ElementCollection& element : mMesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        grpCellsMT.Add(cellContainer.back());
    }

    // create Boundary elements
    InterpolationSimple& boundaryInterpolation = mMesh.CreateInterpolation(InterpolationPoint());

    Eigen::VectorXd coordLeftBoundary = Eigen::VectorXd::Zero(1);
    Eigen::VectorXd coordRightBoundary = Eigen::VectorXd::Ones(1);

    ElementCollectionFem& leftElementCollection =
            mMesh.Elements.Add({{{mMesh.NodeAtCoordinate(coordLeftBoundary)}, boundaryInterpolation}});
    leftElementCollection.AddDofElement(
            dofWV, ElementFem({&mMesh.NodeAtCoordinate(coordLeftBoundary, dofWV)}, boundaryInterpolation));
    leftElementCollection.AddDofElement(
            dofRH, ElementFem({&mMesh.NodeAtCoordinate(coordLeftBoundary, dofRH)}, boundaryInterpolation));


    ElementCollectionFem& rightElementCollection =
            mMesh.Elements.Add({{{mMesh.NodeAtCoordinate(coordRightBoundary)}, boundaryInterpolation}});
    rightElementCollection.AddDofElement(
            dofWV, ElementFem({&mMesh.NodeAtCoordinate(coordRightBoundary, dofWV)}, boundaryInterpolation));
    rightElementCollection.AddDofElement(
            dofRH, ElementFem({&mMesh.NodeAtCoordinate(coordRightBoundary, dofRH)}, boundaryInterpolation));

    cellContainer.push_back(new CellPoint(leftElementCollection, cellId++));
    grpCellsMTBoundary.Add(cellContainer.back());
    cellContainer.push_back(new CellPoint(rightElementCollection, cellId++));
    grpCellsMTBoundary.Add(cellContainer.back());
}


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofVector<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Gradient()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildVector(
            grpCellsMT, {dofRH, dofWV}, std::bind(&MoistureTransport::Gradient, &integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofMatrixSparse<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Stiffness()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildMatrix(
            grpCellsMT, {dofRH, dofWV}, std::bind(&MoistureTransport::Stiffness, &integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofMatrixSparse<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Damping()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildMatrix(
            grpCellsMT, {dofRH, dofWV}, std::bind(&MoistureTransport::Damping, &integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofVector<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::GradientBoundary()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildVector(grpCellsMTBoundary, {dofRH, dofWV},
                                                std::bind(&MoistureTransportBoundary<TDim, TDCw, TDCg, TWVEq>::Gradient,
                                                          mIntegrandMoistureTransportBoundary, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofMatrixSparse<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::StiffnessBoundary()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildMatrix(
            grpCellsMTBoundary, {dofRH, dofWV},
            std::bind(&MoistureTransportBoundary<TDim, TDCw, TDCg, TWVEq>::Stiffness,
                      mIntegrandMoistureTransportBoundary, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
DofVector<double> MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::CreateDofVector()
{
    DofVector<double> d;

    d[dofWV] = Eigen::VectorXd(mMesh.NodesTotal(dofWV).Size() * dofWV.GetNum());
    d[dofRH] = Eigen::VectorXd(mMesh.NodesTotal(dofRH).Size() * dofRH.GetNum());

    return d;
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::GetDofVector(DofVector<double>& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Extract(&dofs, dofs.DofTypes(), instance);
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::ExtractDofs()
{
    DofVector<double> dg{CreateDofVector()}, vg{CreateDofVector()}, ag{CreateDofVector()};
    GetDofVector(dg);
    GetDofVector(vg, 1);
    GetDofVector(ag, 2);
    return std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>(ToEigen(dg, GetDofs()), ToEigen(vg, GetDofs()),
                                                                         ToEigen(ag, GetDofs()));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MergeDofVector(DofVector<double>& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Merge(dofs, dofs.DofTypes(), instance);
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MergeDofs(Eigen::VectorXd d, Eigen::VectorXd v,
                                                                     Eigen::VectorXd a)
{
    DofVector<double> dg{CreateDofVector()}, vg{CreateDofVector()}, ag{CreateDofVector()};
    FromEigen(d, GetDofs(), &dg);
    FromEigen(v, GetDofs(), &vg);
    FromEigen(a, GetDofs(), &ag);
    MergeDofVector(dg);
    MergeDofVector(vg, 1);
    MergeDofVector(ag, 2);
}
