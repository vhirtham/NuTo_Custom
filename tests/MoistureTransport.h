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

// NuToCustom includes
#include "integrands/MoistureTransport.h"

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


    GlobalDofVector Gradient();
    GlobalDofMatrixSparse Stiffness();
    GlobalDofMatrixSparse Damping();


    GlobalDofVector CreateDofVector();
    void GetDofVector(GlobalDofVector& dofs, int instance = 0);
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ExtractDofs();
    void MergeDofVector(GlobalDofVector& dofs, int instance = 0);
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
    IntegrationTypeTensorProduct<1> integrationType = {3, eIntegrationMethod::GAUSS};
    MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq> integrandMoistureTransport;
};


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MoistureTransportTest(double rho_w, double rho_g_sat, double PV)
    : integrandMoistureTransport(dofWV, dofRH, rho_w, rho_g_sat, PV)
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
}


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
GlobalDofVector MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Gradient()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildVector(
            grpCellsMT, {dofRH, dofWV},
            std::bind(&MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Gradient, integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
GlobalDofMatrixSparse MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Stiffness()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildMatrix(grpCellsMT, {dofRH, dofWV},
                                                std::bind(&MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Stiffness,
                                                          integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
GlobalDofMatrixSparse MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::Damping()
{
    CheckDofNumbering();
    return SimpleAssembler(dofInfo).BuildMatrix(
            grpCellsMT, {dofRH, dofWV},
            std::bind(&MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Damping, integrandMoistureTransport, _1, 0.));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
GlobalDofVector MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::CreateDofVector()
{
    GlobalDofVector d;
    int nDofWV = (mMesh.NodesTotal(dofWV)).Size();
    int nDofWV_dep = constraints.GetNumEquations(dofWV);
    int nDofRH = (mMesh.NodesTotal(dofRH)).Size();
    int nDofRH_dep = constraints.GetNumEquations(dofRH);

    d.J[dofWV] = Eigen::VectorXd(nDofWV - nDofWV_dep);
    d.K[dofWV] = Eigen::VectorXd(nDofWV_dep);
    d.J[dofRH] = Eigen::VectorXd(nDofRH - nDofRH_dep);
    d.K[dofRH] = Eigen::VectorXd(nDofRH_dep);

    return d;
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::GetDofVector(GlobalDofVector& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Extract(&dofs, dofs.J.DofTypes(), instance);
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::ExtractDofs()
{
    GlobalDofVector dg{CreateDofVector()}, vg{CreateDofVector()}, ag{CreateDofVector()};
    GetDofVector(dg);
    GetDofVector(vg, 1);
    GetDofVector(ag, 2);
    return std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>(
            ToEigen(dg.J, GetDofs()), ToEigen(vg.J, GetDofs()), ToEigen(ag.J, GetDofs()));
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MergeDofVector(GlobalDofVector& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Merge(dofs, dofs.J.DofTypes(), instance);
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MergeDofs(Eigen::VectorXd d, Eigen::VectorXd v,
                                                                     Eigen::VectorXd a)
{
    GlobalDofVector dg{CreateDofVector()}, vg{CreateDofVector()}, ag{CreateDofVector()};
    FromEigen(d, GetDofs(), &dg.J);
    FromEigen(v, GetDofs(), &vg.J);
    FromEigen(a, GetDofs(), &ag.J);
    MergeDofVector(dg);
    MergeDofVector(vg, 1);
    MergeDofVector(ag, 2);
}
