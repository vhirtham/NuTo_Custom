#pragma once

#include <cmath>
#include <memory>
namespace NuTo
{

class MTCoefficientInterface
{
public:
    virtual std::unique_ptr<MTCoefficientInterface> Clone() const = 0;
    virtual double value(double WV, double RH, double WV_dt, double RH_dt) const = 0;
    virtual double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt) const = 0;
    virtual double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt) const = 0;
    virtual double d_WaterVolumeFraction_dt(double WV, double RH, double WV_dt, double RH_dt) const = 0;
    virtual double d_RelativeHumidity_dt(double WV, double RH, double WV_dt, double RH_dt) const = 0;
};


class MTCoefficientConstant : public MTCoefficientInterface
{
    double mC0;

public:
    MTCoefficientConstant(double c0)
        : mC0{c0}
    {
    }

    virtual std::unique_ptr<MTCoefficientInterface> Clone() const final
    {
        return std::make_unique<MTCoefficientConstant>(*this);
    }

    virtual double value(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC0;
    }

    virtual double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return 0.;
    }

    virtual double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return 0.;
    }

    virtual double d_WaterVolumeFraction_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return 0.;
    }

    virtual double d_RelativeHumidity_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return 0.;
    }
};

class MTCoefficientLinear : public MTCoefficientInterface
{
    double mC1WV, mC1RH, mC1WV_dt, mC1RH_dt, mC0;

public:
    MTCoefficientLinear(double c1WV, double c1RH, double c1WV_dt, double c1RH_dt, double c0)
        : mC1WV{c1WV}
        , mC1RH{c1RH}
        , mC1WV_dt{c1WV_dt}
        , mC1RH_dt{c1RH_dt}
        , mC0{c0}
    {
    }

    virtual std::unique_ptr<MTCoefficientInterface> Clone() const final
    {
        return std::make_unique<MTCoefficientLinear>(*this);
    }

    virtual double value(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC1WV * WV + mC1RH * RH + mC1WV_dt * WV_dt + mC1RH_dt * RH_dt + mC0;
    }

    virtual double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC1WV;
    }

    virtual double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC1RH;
    }

    virtual double d_WaterVolumeFraction_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC1WV_dt;
    }

    virtual double d_RelativeHumidity_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return mC1RH_dt;
    }
};


class MTCoefficientQuadratic : public MTCoefficientInterface
{
    std::array<double, 2> mC1WV, mC1RH, mC1WV_dt, mC1RH_dt;
    double mC0;

public:
    MTCoefficientQuadratic(std::array<double, 2> c1WV, std::array<double, 2> c1RH, std::array<double, 2> c1WV_dt,
                           std::array<double, 2> c1RH_dt, double c0)
        : mC1WV{c1WV}
        , mC1RH{c1RH}
        , mC1WV_dt{c1WV_dt}
        , mC1RH_dt{c1RH_dt}
        , mC0{c0}
    {
    }

    virtual std::unique_ptr<MTCoefficientInterface> Clone() const final
    {
        return std::make_unique<MTCoefficientQuadratic>(*this);
    }

    virtual double value(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return dofValue(mC1WV, WV) + dofValue(mC1RH, RH) + dofValue(mC1WV_dt, WV_dt) + dofValue(mC1RH_dt, RH_dt) + mC0;
    }

    virtual double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return dofGradient(mC1WV, WV);
    }

    virtual double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return dofGradient(mC1RH, RH);
    }

    virtual double d_WaterVolumeFraction_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return dofGradient(mC1WV_dt, WV_dt);
    }

    virtual double d_RelativeHumidity_dt(double WV, double RH, double WV_dt, double RH_dt) const final
    {
        return dofGradient(mC1RH_dt, RH_dt);
    }

private:
    inline double dofValue(std::array<double, 2> Coeffs, double dof) const
    {
        return Coeffs[1] * dof * dof + Coeffs[0] * dof;
    }

    inline double dofGradient(std::array<double, 2> Coeffs, double dof) const
    {
        return 2 * Coeffs[1] * dof + Coeffs[0];
    }
};
}
