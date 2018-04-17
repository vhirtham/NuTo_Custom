#include "structures/MultiPhysicsStructure.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/tools/NodalValueMerger.h"

using namespace std::placeholders;

namespace NuTo
{

DofVector<double> MultiPhysicsStructure::GradientMoistureTransport()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildVector(
            mGrpCellsMT, {mDofRH, mDofWV},
            std::bind(&MoistureTransportIntegrand::Gradient, mMoistureTransport, _1, 0.));
}

DofVector<double> MultiPhysicsStructure::GradientMoistureTransportBoundary()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildVector(
            mGrpCellsMTBoundary, {mDofRH, mDofWV},
            std::bind(&MoistureTransportBoundaryIntegrand::Gradient, mMoistureTransportBoundary, _1, 0.));
}

DofMatrixSparse<double> MultiPhysicsStructure::StiffnessMoistureTransport()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildMatrix(
            mGrpCellsMT, {mDofRH, mDofWV},
            std::bind(&MoistureTransportIntegrand::Stiffness, mMoistureTransport, _1, 0.));
}

DofMatrixSparse<double> MultiPhysicsStructure::StiffnessMoistureTransportBoundary()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildMatrix(
            mGrpCellsMTBoundary, {mDofRH, mDofWV},
            std::bind(&MoistureTransportBoundaryIntegrand::Stiffness, mMoistureTransportBoundary, _1, 0.));
}

DofMatrixSparse<double> MultiPhysicsStructure::DampingMoistureTransport()
{
    CheckDofNumbering();
    return SimpleAssembler(mDofInfo).BuildMatrix(
            mGrpCellsMT, {mDofRH, mDofWV}, std::bind(&MoistureTransportIntegrand::Damping, mMoistureTransport, _1, 0.));
}

void MultiPhysicsStructure::CheckDofNumbering()
{
    if (!(mDofInfo.numDependentDofs.Has(mDofDisp) && mDofInfo.numDependentDofs.Has(mDofWV) &&
          mDofInfo.numDependentDofs.Has(mDofRH)))
        throw Exception(__PRETTY_FUNCTION__, "Please do a dof numbering before trying to assemble something!");
}

DofVector<double> MultiPhysicsStructure::CreateGlobalDofVector(std::vector<DofType> dofs)
{
    DofVector<double> d;
    for (const auto& dof : dofs)
    {
        d[dof] = Eigen::VectorXd(mMesh.NodesTotal(dof).Size() * dof.GetNum());
    }
    return d;
}

void MultiPhysicsStructure::GetDofVector(DofVector<double>& dofVector, std::vector<DofType> dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Extract(&dofVector, dofs, instance);
}

void MultiPhysicsStructure::MergeDofs(Eigen::VectorXd values, std::vector<DofType> dofs, int instance)
{
    DofVector<double> gdv{CreateGlobalDofVector(dofs)};
    FromEigen(values, dofs, &gdv);
    MergeDofVector(gdv, instance);
}

void MultiPhysicsStructure::MergeDofVector(DofVector<double>& dofs, int instance)
{
    CheckDofNumbering();
    NodalValueMerger Merger(&mMesh);
    Merger.Merge(dofs, dofs.DofTypes(), instance);
}

void MultiPhysicsStructure::CreateUnitMesh(int numX, int numY, const InterpolationSimple& interpolationDispArg,
                                           const InterpolationSimple& interpolationWVArg,
                                           const InterpolationSimple& interpolationRHArg,
                                           const InterpolationSimple& interpolationMTBoundaryArg)
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

    InterpolationSimple& boundaryInterpolation = mMesh.CreateInterpolation(interpolationMTBoundaryArg);
    AddBoundaryElementsAtAxis(cellId, boundaryInterpolation, eDirection::X);
    AddBoundaryElementsAtAxis(cellId, boundaryInterpolation, eDirection::X, 1.);
    AddBoundaryElementsAtAxis(cellId, boundaryInterpolation, eDirection::Y);
    AddBoundaryElementsAtAxis(cellId, boundaryInterpolation, eDirection::Y, 1.);
}

void MultiPhysicsStructure::AddBoundaryElementsAtAxis(int& cellId, const InterpolationSimple& interpolationMTBoundary,
                                                      eDirection direction, double axisOffset)
{
    // create Boundary elements


    auto grpBoundaryNodeCoords = mMesh.NodesAtAxis(direction, axisOffset);

    std::vector<Eigen::VectorXd> NodeCoords;
    for (auto node : grpBoundaryNodeCoords)
    {
        NodeCoords.push_back(node.GetValues());
    }

    int directionIndex = 0;
    if (direction == eDirection::X)
        directionIndex = 1;


    std::sort(NodeCoords.begin(), NodeCoords.end(), [=](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
        return a[directionIndex] < b[directionIndex];
    });


    for (unsigned int i = 0; i < NodeCoords.size() - 2; i += 2)
    {


        ElementCollectionFem& lowerElementCollection =
                mMesh.Elements.Add({{{mMesh.NodeAtCoordinate(NodeCoords[i]), mMesh.NodeAtCoordinate(NodeCoords[i + 1]),
                                      mMesh.NodeAtCoordinate(NodeCoords[i + 2])},
                                     interpolationMTBoundary}});
        lowerElementCollection.AddDofElement(mDofWV, ElementFem({&mMesh.NodeAtCoordinate(NodeCoords[i], mDofWV),
                                                                 &mMesh.NodeAtCoordinate(NodeCoords[i + 1], mDofWV),
                                                                 &mMesh.NodeAtCoordinate(NodeCoords[i + 2], mDofWV)},
                                                                interpolationMTBoundary));
        lowerElementCollection.AddDofElement(mDofRH, ElementFem({&mMesh.NodeAtCoordinate(NodeCoords[i], mDofRH),
                                                                 &mMesh.NodeAtCoordinate(NodeCoords[i + 1], mDofRH),
                                                                 &mMesh.NodeAtCoordinate(NodeCoords[i + 2], mDofRH)},
                                                                interpolationMTBoundary));

        mCellContainer.push_back(new Cell(lowerElementCollection, mIntegrationTypeBoundary, cellId++));
        mGrpCellsMTBoundary.Add(mCellContainer.back());
    }
}

Eigen::VectorXd MultiPhysicsStructure::ExtractDofs(std::vector<DofType> dofs, int instance)
{
    DofVector<double> gdv{CreateGlobalDofVector(dofs)};
    GetDofVector(gdv, dofs, instance);
    return ToEigen(gdv, dofs);
}
}
