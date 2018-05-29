#include "MultiPhysicsStructureNew.h"

#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/dofs/DofNumbering.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/tools/NodalValueMerger.h"

#include <cassert>

using namespace std::placeholders;

namespace NuTo
{

MultiPhysicsStructure::MultiPhysicsStructure(MeshFem& mesh)
    : mMesh{mesh}
{
}

const DofType& NuTo::MultiPhysicsStructure::AddDofType(std::string name, int numValues)
{
    mDofTypes.emplace_back(name, numValues);
    return mDofTypes.back();
}

Group<CellInterface> MultiPhysicsStructure::CreateCells(const NuTo::Group<ElementCollectionFem>& elementCollection,
                                                        const IntegrationTypeBase& integrationType)
{
    Group<CellInterface> cellGroup;
    for (const auto& element : elementCollection)
    {
        mCells.emplace_back(new Cell(element, integrationType, GetCellId()));
        cellGroup.Add(*mCells.back());
    }
    return cellGroup;
}

Group<CellInterface> MultiPhysicsStructure::GetAllCells() const
{
    Group<CellInterface> cellGroup;
    for (const auto& cell : mCells)
        cellGroup.Add(*cell);
    return cellGroup;
}

const DofInfo& MultiPhysicsStructure::GetDofInfo() const
{
    return mDofInfo;
}


void MultiPhysicsStructure::RenumberDofs()
{
    for (const auto& dofType : mDofTypes)
        mDofInfo.Merge(dofType, DofNumbering::Build(mMesh.NodesTotal(dofType), dofType, mConstraints));
}

void MultiPhysicsStructure::SetNodalValues(const DofType& dofType, std::vector<double> values, int instance)
{
    SetNodalValues(dofType, [&values, instance](NodeSimple& node) {
        for (unsigned int i = 0; i < values.size(); ++i)
            node.SetValue(i, values[i], instance);
    });
}

void MultiPhysicsStructure::SetNodalValues(const DofType& dofType, std::function<void(NodeSimple&)> valueFunction)
{
    for (NodeSimple& node : mMesh.NodesTotal(dofType))
        valueFunction(node);
}

void MultiPhysicsStructure::SetNumTimeDerivatives(int numTimeDerivatives)
{
    assert(numTimeDerivatives > -1);
    mNumTimeDerivatives = numTimeDerivatives;
    for (const auto& dofType : mDofTypes)
        mMesh.AllocateDofInstances(dofType, mNumTimeDerivatives + 1);
}

Eigen::VectorXd MultiPhysicsStructure::AssembleResidual(
        std::vector<DofType> dofTypes,
        std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs)
{
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(0);

    for (const auto& cellGroupIntegrandPair : cellGroupIntegrandPairs)
    {
        assert(cellGroupIntegrandPair.first != nullptr && cellGroupIntegrandPair.second != nullptr);

        auto dofVector = SimpleAssembler(mDofInfo).BuildVector(
                *cellGroupIntegrandPair.first, dofTypes,
                std::bind(&IntegrandWrapperInterface<double>::Residual, cellGroupIntegrandPair.second, _1, 0.));
        if (residual.rows() == 0)
            residual = ToEigen(dofVector, dofTypes);
        else
            residual += ToEigen(dofVector, dofTypes);
    }
    return residual;
}

Eigen::SparseMatrix<double> MultiPhysicsStructure::AssembleStiffness(
        std::vector<DofType> dofTypes,
        std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs)
{
    return AssembleMatrix(dofTypes, cellGroupIntegrandPairs, &IntegrandWrapperInterface<double>::Stiffness);
}

Eigen::SparseMatrix<double> MultiPhysicsStructure::AssembleDamping(
        std::vector<DofType> dofTypes,
        std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs)
{
    return AssembleMatrix(dofTypes, cellGroupIntegrandPairs, &IntegrandWrapperInterface<double>::Damping);
}

Eigen::SparseMatrix<double> MultiPhysicsStructure::AssembleMass(
        std::vector<DofType> dofTypes,
        std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs)
{
    return AssembleMatrix(dofTypes, cellGroupIntegrandPairs, &IntegrandWrapperInterface<double>::Mass);
}

int MultiPhysicsStructure::GetCellId()
{
    static int id = 0;
    return id++;
}

Eigen::SparseMatrix<double> MultiPhysicsStructure::AssembleMatrix(
        std::vector<DofType> dofTypes,
        std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs,
        DofMatrix<double> (IntegrandWrapperInterface<double>::*matrixFunction)(const CellIpData&, double))
{
    Eigen::SparseMatrix<double> matrix;
    assert(matrix.rows() == 0 && matrix.cols() == 0);

    for (const auto& cellGroupIntegrandPair : cellGroupIntegrandPairs)
    {
        assert(cellGroupIntegrandPair.first != nullptr && cellGroupIntegrandPair.second != nullptr);
        auto dofMatrix =
                SimpleAssembler(mDofInfo).BuildMatrix(*cellGroupIntegrandPair.first, dofTypes,
                                                      std::bind(matrixFunction, cellGroupIntegrandPair.second, _1, 0.));
        if (matrix.rows() == 0)
            matrix = ToEigen(dofMatrix, dofTypes);
        else
            matrix += ToEigen(dofMatrix, dofTypes);
    }
    return matrix;
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

Eigen::VectorXd MultiPhysicsStructure::ExtractDofs(std::vector<DofType> dofs, int instance)
{
    if (!IsDofNumberingDone())
        throw Exception(__PRETTY_FUNCTION__, "Please do a DOF numbering first");
    auto dofVector = CreateGlobalDofVector(dofs);
    NodalValueMerger Merger(&mMesh);
    Merger.Extract(&dofVector, dofs, instance);
    return ToEigen(dofVector, dofs);
}

void MultiPhysicsStructure::MergeDofs(Eigen::VectorXd values, std::vector<DofType> dofs, int instance)
{
    if (!IsDofNumberingDone())
        throw Exception(__PRETTY_FUNCTION__, "Please do a DOF numbering first");
    auto dofVector = CreateGlobalDofVector(dofs);
    FromEigen(values, dofs, &dofVector);
    NodalValueMerger Merger(&mMesh);
    Merger.Merge(dofVector, dofs, instance);
}

bool MultiPhysicsStructure::IsDofNumberingDone() const
{
    for (const auto& dofType : mDofTypes)
        if (!mDofInfo.numDependentDofs.Has(dofType))
            return false;
    return true;
}
}
