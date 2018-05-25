#include "MultiPhysicsStructureNew.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/dofs/DofNumbering.h"

#include <cassert>

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

int MultiPhysicsStructure::GetCellId()
{
    static int id = 0;
    return id++;
}
}
