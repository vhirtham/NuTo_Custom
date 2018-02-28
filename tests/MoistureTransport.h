#pragma once

// NuTo includes
#include "mechanics/cell/SimpleAssembler.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"

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
    MoistureTransportTest(double rho_w = 1.0, double rho_g_sat = 1.0, double PV = 0.2);

    void CreateUnitMesh(int numElements, const InterpolationSimple& interpolationWVArg,
                        const InterpolationSimple& interpolationRHArg);


    GlobalDofVector Gradient();
    GlobalDofMatrixSparse Stiffness();
    GlobalDofMatrixSparse Damping();


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
    MeshFem mesh;
    boost::ptr_vector<CellInterface> cellContainer;
    Group<CellInterface> grpCellsMT;
    IntegrationTypeTensorProduct<1> integrationType = {3, eIntegrationMethod::GAUSS};
    MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq> integrandMoistureTransport;
};


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::MoistureTransportTest(double rho_w, double rho_g_sat, double PV)
    : integrandMoistureTransport(dofRH, dofWV, rho_w, rho_g_sat, PV)
{
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void MoistureTransportTest<TDim, TDCw, TDCg, TMeC, TWVEq>::CreateUnitMesh(int numElements,
                                                                          const InterpolationSimple& interpolationWVArg,
                                                                          const InterpolationSimple& interpolationRHArg)
{
    // mesh creation
    mesh = UnitMeshFem::CreateLines(numElements);
    const auto& interpolationWV = mesh.CreateInterpolation(interpolationWVArg);
    const auto& interpolationRH = mesh.CreateInterpolation(interpolationRHArg);
    AddDofInterpolation(&mesh, dofWV, interpolationWV);
    AddDofInterpolation(&mesh, dofRH, interpolationRH);

    // dof numbering
    dofInfo = DofNumbering::Build(mesh.NodesTotal(dofRH), dofRH, constraints);
    dofInfo.Merge(dofWV, DofNumbering::Build(mesh.NodesTotal(dofWV), dofWV, constraints));

    // set nodal values to 100% RH and equilibium state
    for (NodeSimple& node : mesh.NodesTotal(dofRH))
        node.SetValue(0, 1.0);
    for (NodeSimple& node : mesh.NodesTotal(dofWV))
        node.SetValue(0, TWVEq::value(0., 1., 0., 0.));


    // create cells
    int cellId = 0;
    for (ElementCollection& element : mesh.Elements)
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
