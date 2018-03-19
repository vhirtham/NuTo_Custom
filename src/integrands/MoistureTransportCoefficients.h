#pragma once

#include <cmath>

namespace NuTo
{


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
