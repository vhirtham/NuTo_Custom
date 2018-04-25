#pragma once
#include "MoistureTransport.h"


#define EXTRACTSHAREDVARS                                                                                              \
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofTypeWV); /*discretized water volume fraction*/           \
    const Eigen::MatrixXd Nw = cellIpData.N(mDofTypeWV); /*shape function of the water phase*/                         \
    const Eigen::MatrixXd Bw = cellIpData.B(mDofTypeWV, Nabla::Gradient()); /* deriv. shape func. of the water phase*/ \
    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofTypeRH); /*discretized relative humidity*/               \
    const Eigen::MatrixXd Ng = cellIpData.N(mDofTypeRH); /*shape function of the gas phase*/                           \
    const Eigen::MatrixXd Bg = cellIpData.B(mDofTypeRH, Nabla::Gradient()); /* deriv. shape func. of the gas phase*/   \
    const double wv = (Nw * WV)[0]; /*scalar water volume fraction*/                                                   \
    const double wv_dt = (Nw * cellIpData.NodeValueVector(mDofTypeWV,1))[0]; /*scalar water volume fraction vel.*/    \
    const double rh = (Ng * RH)[0]; /*scalar relative humidity*/                                                       \
    const double rh_dt = (Ng * cellIpData.NodeValueVector(mDofTypeRH,1))[0]; /*scalar relative humidity vel.*/        \
    const double mec = TMeC::value(wv, rh, wv_dt, rh_dt); /*mass exchange coefficient*/                                \
    const double wveq = TWVEq::value(wv, rh, wv_dt, rh_dt); /*equilibrium water volume fraction*/

#define EXTRACTGRADIENTVARS                                                                                            \
    EXTRACTSHAREDVARS                                                                                                  \
    const double dcw = TDCw::value(wv, rh, wv_dt, rh_dt); /*diffusion coefficient of the water phase*/                 \
    const double dcg = TDCg::value(wv, rh, wv_dt, rh_dt); /*diffusion coefficient of the gas phase*/

#define EXTRACTSTIFFNESSVARS                                                                                           \
    EXTRACTGRADIENTVARS                                                                                                \
    const Eigen::MatrixXd Iw = Eigen::MatrixXd::Identity(WV.rows(), WV.rows()); /*identity matrix of the water phase*/ \
    const Eigen::MatrixXd Ig = Eigen::MatrixXd::Identity(RH.rows(), RH.rows()); /*identity matrix of the gas phase*/   \
    /*derivatives with the naming pattern VAR1_dVAR2*/                                                                 \
    const double dcw_dwv = TDCw::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double dcw_drh = TDCw::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double dcg_dwv = TDCg::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double dcg_drh = TDCg::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double mec_dwv = TMeC::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double mec_drh = TMeC::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double wveq_dwv = TWVEq::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                        \
    const double wveq_drh = TWVEq::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);

#define EXTRACTDAMPINGVARS                                                                                             \
    EXTRACTSHAREDVARS                                                                                                  \
    const Eigen::MatrixXd GWV = Bw * WV; /*water volume graction gradient*/                                            \
    const Eigen::MatrixXd GRH = Bg * RH; /*relative humidity gradient*/                                                \
    /*derivatives with the naming pattern VAR1_dVAR2*/                                                                 \
    const double dcw_dwv_dt = TDCw::d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double dcw_drh_dt = TDCw::d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double dcg_dwv_dt = TDCg::d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double dcg_drh_dt = TDCg::d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double mec_dwv_dt = TMeC::d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double mec_drh_dt = TMeC::d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double wveq_dwv_dt = TWVEq::d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                  \
    const double wveq_drh_dt = TWVEq::d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);


template <typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::Integrands::MoistureTransport<TDCw, TDCg, TMeC, TWVEq>::MoistureTransport(DofType dofTypeWV, DofType dofTypeRH, double rho_w,
        double rho_g_sat, double PV)
    : mDofTypeWV(dofTypeWV)
    , mDofTypeRH(dofTypeRH)
    , mRho_w(rho_w)
    , mRho_g_sat(rho_g_sat)
    , mPV(PV)
{
}

template <typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofVector<double>
NuTo::Integrands::MoistureTransport<TDCw, TDCg, TMeC, TWVEq>::Gradient(const NuTo::CellIpData& cellIpData,
                                                                             double deltaT)
{
    DofVector<double> gradient;

    // variable meanings can be found in macro definition at the beginning of the file
    EXTRACTGRADIENTVARS
    CheckValuesValid(wv);

    // superpositions
    const double me = (wveq - wv) * mec;

    // clang-format off
    gradient[mDofTypeWV] = -dcw * Bw.transpose() * Bw * WV
                         + (me - mRho_w * wv_dt) * Nw.transpose();

    gradient[mDofTypeRH] = -dcg * Bg.transpose() * Bg * RH
                         + (mRho_g_sat * ((wv - mPV) * rh_dt + rh * wv_dt) - me) * Ng.transpose();
    // clang-format on


    return gradient;
}

template <typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofMatrix<double>
NuTo::Integrands::MoistureTransport<TDCw, TDCg, TMeC, TWVEq>::Stiffness(const NuTo::CellIpData& cellIpData,
                                                                              double deltaT)
{
    NuTo::DofMatrix<double> stiffness;

    // variable meanings can be found in macro definition at the beginning of the file
    EXTRACTSTIFFNESSVARS
    CheckValuesValid(wv);

    // superpositions
    double supA = mec_dwv * (wveq - wv) + mec * wveq_dwv - mec;
    double supB = mec_drh * (wveq - wv) + mec * wveq_drh;

    // clang-format off
    stiffness(mDofTypeWV, mDofTypeWV) = -Bw.transpose() * Bw * (dcw * Iw + dcw_dwv * WV * Nw)
                                      +  Nw.transpose() * supA * Nw;

    stiffness(mDofTypeWV, mDofTypeRH) = -dcw_drh * Bw.transpose() * Bw * WV * Ng
                                      +  Nw.transpose() * supB * Ng;

    stiffness(mDofTypeRH, mDofTypeWV) = -dcg_dwv * Bg.transpose() * Bg * RH * Nw
                                      +  Ng.transpose() * (-supA + mRho_g_sat * rh_dt) * Nw;

    stiffness(mDofTypeRH, mDofTypeRH) = -Bg.transpose() * Bg * (dcg * Ig + dcg_drh * RH * Ng)
                                      +  Ng.transpose() * (-supB + mRho_g_sat * wv_dt) * Ng;
    // clang-format on


    return stiffness;
}

template <typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofMatrix<double>
NuTo::Integrands::MoistureTransport<TDCw, TDCg, TMeC, TWVEq>::Damping(const NuTo::CellIpData& cellIpData,
                                                                            double deltaT)
{
    NuTo::DofMatrix<double> damping;

    // variable meanings can be found in macro definition at the beginning of the file
    EXTRACTDAMPINGVARS
    CheckValuesValid(wv);
    // superpositions
    const double supA = mec_dwv_dt * (wveq - wv) + mec * wveq_dwv_dt;
    const double supB = mec_drh_dt * (wveq - wv) + mec * wveq_drh_dt;

    // clang-format off
    damping(mDofTypeWV, mDofTypeWV) = -dcw_dwv_dt * Bw.transpose() * GWV * Nw
                                    + Nw.transpose() * (supA - mRho_w) * Nw;

    damping(mDofTypeWV, mDofTypeRH) = -dcw_drh_dt * Bw.transpose() * GWV * Ng
                                    + Nw.transpose() * supB * Ng;

    damping(mDofTypeRH, mDofTypeWV) = -dcg_dwv_dt * Bg.transpose() * GRH * Nw
                                    + Ng.transpose() * (mRho_g_sat * rh - supA) * Nw;

    damping(mDofTypeRH, mDofTypeRH) = -dcg_drh_dt * Bg.transpose() * GRH * Ng
                                    + Ng.transpose() * (mRho_g_sat * (wv - mPV) - supB) * Ng;
    // clang-format on


    return damping;
}

template < typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
void NuTo::Integrands::MoistureTransport<TDCw, TDCg, TMeC, TWVEq>::CheckValuesValid(
        double wv)
{
    if (std::abs(mPV - wv) < 10e-9)
        throw Exception(
                __PRETTY_FUNCTION__,
                "Poresvolume is completly occupied by water. Gas transport equation will produce undefined behaviour!");
}


#undef EXTRACTDAMPINGVARS
#undef EXTRACTSTIFFNESSVARS
#undef EXTRACTGRADIENTVARS
#undef EXTRACTSHAREDVARS
