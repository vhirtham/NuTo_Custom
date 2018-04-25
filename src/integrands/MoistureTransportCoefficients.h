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


template <int TC1WV, int TC1RH, int TC1WV_dt, int TC1RH_dt, int TC0, int TDiv = 1>
class MTCLinear
{
    constexpr static double mC0 = static_cast<double>(TC0) / static_cast<double>(TDiv);
    constexpr static double mC1WV = static_cast<double>(TC1WV) / static_cast<double>(TDiv);
    constexpr static double mC1RH = static_cast<double>(TC1RH) / static_cast<double>(TDiv);
    constexpr static double mC1WV_dt = static_cast<double>(TC1WV_dt) / static_cast<double>(TDiv);
    constexpr static double mC1RH_dt = static_cast<double>(TC1RH_dt) / static_cast<double>(TDiv);

public:
    constexpr static double value(double WV, double RH, double WV_dt, double RH_dt)
    {
        return mC1WV * WV + mC1RH * RH + mC1WV_dt * WV_dt + mC1RH_dt * RH_dt + mC0;
    }
    constexpr static double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt)
    {
        return mC1WV;
    }
    constexpr static double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt)
    {
        return mC1RH;
    }
    constexpr static double d_WaterVolumeFraction_dt(double WV, double RH, double WV_dt, double RH_dt)
    {
        return mC1WV_dt;
    }
    constexpr static double d_RelativeHumidity_dt(double WV, double RH, double WV_dt, double RH_dt)
    {
        return mC1RH_dt;
    }
};

// Linear functions depending on only one dof
template <int TC1, int TC0, int TDiv = 1>
using MTCLinearWV = MTCLinear<TC1, 0, 0, 0, TC0, TDiv>;
template <int TC1, int TC0, int TDiv = 1>
using MTCLinearRH = MTCLinear<0, TC1, 0, 0, TC0, TDiv>;
template <int TC1, int TC0, int TDiv = 1>
using MTCLinearWV_dt = MTCLinear<0, 0, TC1, 0, TC0, TDiv>;
template <int TC1, int TC0, int TDiv = 1>
using MTCLinearRH_dt = MTCLinear<0, 0, 0, TC1, TC0, TDiv>;

// Constant functions
template <int TC0, int TDiv = 1>
using MTCConst = MTCLinear<0, 0, 0, 0, TC0, TDiv>;
using MTCZero = MTCLinear<0, 0, 0, 0, 0, 0>;
}
