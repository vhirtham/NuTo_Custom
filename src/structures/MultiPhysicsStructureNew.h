#pragma once

#include "nuto/base/Group.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

#include "tools/IntegrandWrapper.h"

#include <list>
#include <functional>
#include <string>
#include <vector>

namespace NuTo
{

class MeshFem;

class MultiPhysicsStructure
{
    MeshFem& mMesh;
    std::list<DofType> mDofTypes;
    DofInfo mDofInfo;
    Constraint::Constraints mConstraints;
    std::vector<std::unique_ptr<IntegrationTypeBase>> mIntegrationTypes;
    int mNumTimeDerivatives = 0;

    // cells
    std::vector<std::unique_ptr<CellInterface>> mCells;

    // integrands
    std::vector<std::unique_ptr<IntegrandWrapperInterface<double>>> mIntegrands;

public:
    MultiPhysicsStructure(MeshFem& mesh);

    const DofType& AddDofType(std::string name, int numValues);

    template <class TIntegrand>
    TIntegrand& AddIntegrand(TIntegrand&& integrand)
    {
        mIntegrands.emplace_back(std::make_unique<TIntegrand>(std::move(integrand)));
        return static_cast<TIntegrand&>(*mIntegrands.back());
    }

    template <typename TIntegrationType>
    const TIntegrationType& AddIntegrationType(const TIntegrationType& integrationType)
    {
        mIntegrationTypes.push_back(std::make_unique<TIntegrationType>(integrationType));
        return static_cast<const TIntegrationType&>(*mIntegrationTypes.back());
    }

    Group<CellInterface> CreateCells(const NuTo::Group<ElementCollectionFem>& elementCollection,
                                     const IntegrationTypeBase& integrationType);

    Group<CellInterface> GetAllCells() const;

    const DofInfo& GetDofInfo() const;
    Constraint::Constraints& GetConstraints();
    MeshFem& GetMesh();

    void RenumberDofs();

    void SetNodalValues(const DofType& dofType, std::vector<double> values, int instance = 0);
    void SetNodalValues(const DofType& dofType, std::function<void(NodeSimple&)> valueFunction);
    void SetNumTimeDerivatives(int numTimeDerivatives);

    Eigen::VectorXd AssembleResidual(
            std::vector<DofType> dofTypes,
            std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs);

    Eigen::SparseMatrix<double> AssembleStiffness(
            std::vector<DofType> dofTypes,
            std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs);

    Eigen::SparseMatrix<double> AssembleDamping(
            std::vector<DofType> dofTypes,
            std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs);

    Eigen::SparseMatrix<double> AssembleMass(
            std::vector<DofType> dofTypes,
            std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs);

    Eigen::VectorXd ExtractDofs(std::vector<DofType> dofs, int instance);
    void MergeDofs(Eigen::VectorXd values, std::vector<DofType> dofs, int instance);

private:
    static int GetCellId();

    Eigen::SparseMatrix<double> AssembleMatrix(
            std::vector<DofType> dofTypes,
            std::vector<std::pair<Group<CellInterface>*, IntegrandWrapperInterface<double>*>> cellGroupIntegrandPairs,
            DofMatrix<double> (IntegrandWrapperInterface<double>::*matrixFunction)(const CellIpData&, double));

    DofVector<double> CreateGlobalDofVector(std::vector<DofType> dofs);


    bool IsDofNumberingDone() const;
};
}
