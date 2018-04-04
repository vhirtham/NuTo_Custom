#pragma once

namespace NuTo
{
class InterpolationPoint : public InterpolationSimple
{
    Pyramid mShape;

public:
    virtual std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationPoint>(*this);
    }

    virtual ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return Eigen::VectorXd::Ones(1);
    }

    virtual DerivativeShapeFunctionsNatural
    GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
    }

    virtual NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return Eigen::VectorXd::Ones(1);
    }

    virtual int GetNumNodes() const override
    {
        return 1;
    }

    virtual const Shape& GetShape() const override
    {
        return mShape;
    }
};
}
