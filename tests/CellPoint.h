#pragma once


#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/Jacobian.h"

namespace NuTo
{
class CellPoint : public CellInterface
{
public:
    CellPoint(const ElementCollection& elements, const int id)
        : mElements(elements)
        , mId(id)
        , mShape(elements.GetShape())
    {
    }

    virtual ~CellPoint() = default;

    virtual double Integrate(ScalarFunction) override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    int Id()
    {
        return mId;
    }

    virtual DofVector<double> Integrate(VectorFunction f) override
    {
        DofVector<double> result;
        CellData cellData(mElements, Id());

        Jacobian jacobian{Eigen::VectorXd(), Eigen::MatrixXd()};
        CellIpData cellipData(cellData, jacobian, Eigen::VectorXd(), 0);
        result += f(cellipData);

        return result;
    }
    virtual DofMatrix<double> Integrate(MatrixFunction f) override
    {
        DofMatrix<double> result;
        CellData cellData(mElements, Id());
        Jacobian jacobian{Eigen::VectorXd(), Eigen::MatrixXd()};
        CellIpData cellipData(cellData, jacobian, Eigen::VectorXd(), 0);
        result += f(cellipData);

        return result;
    }
    virtual void Apply(VoidFunction) override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    virtual std::vector<Eigen::VectorXd> Eval(EvalFunction f) const override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    virtual Eigen::VectorXi DofNumbering(DofType dof) override
    {
        return mElements.DofElement(dof).GetDofNumbering();
    }

    //! Coordinate interpolation
    virtual Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords) const override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }
    //! Dof interpolation
    virtual Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords, DofType dof) const override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    virtual const Shape& GetShape() const override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

private:
    const ElementCollection& mElements;
    const int mId;
    const Shape& mShape;
};
}
