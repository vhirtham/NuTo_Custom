#include "MoistureTransport.h"

#include "integrands/MoistureTransportCoefficients.h"

#define EXTRACTSHAREDVARS                                                                                              \
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofTypeWV); /*discretized water volume fraction*/           \
    const Eigen::MatrixXd Nw = cellIpData.N(mDofTypeWV); /*shape function of the water phase*/                         \
    const Eigen::MatrixXd Bw = cellIpData.B(mDofTypeWV, Nabla::Gradient()); /* deriv. shape func. of the water phase*/ \
    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofTypeRH); /*discretized relative humidity*/               \
    const Eigen::MatrixXd Ng = cellIpData.N(mDofTypeRH); /*shape function of the gas phase*/                           \
    const Eigen::MatrixXd Bg = cellIpData.B(mDofTypeRH, Nabla::Gradient()); /* deriv. shape func. of the gas phase*/   \
    const double wv = (Nw * WV)[0]; /*scalar water volume fraction*/                                                   \
    const double wv_dt = (Nw * cellIpData.NodeValueVector(mDofTypeWV, 1))[0]; /*scalar water volume fraction vel.*/    \
    const double rh = (Ng * RH)[0]; /*scalar relative humidity*/                                                       \
    const double rh_dt = (Ng * cellIpData.NodeValueVector(mDofTypeRH, 1))[0]; /*scalar relative humidity vel.*/        \
    const double mec = mMeC->value(wv, rh, wv_dt, rh_dt); /*mass exchange coefficient*/                                \
    const double wveq = mWVEq->value(wv, rh, wv_dt, rh_dt); /*equilibrium water volume fraction*/

#define EXTRACTGRADIENTVARS                                                                                            \
    EXTRACTSHAREDVARS                                                                                                  \
    const double dcw = mDCw->value(wv, rh, wv_dt, rh_dt); /*diffusion coefficient of the water phase*/                 \
    const double dcg = mDCg->value(wv, rh, wv_dt, rh_dt); /*diffusion coefficient of the gas phase*/

#define EXTRACTSTIFFNESSVARS                                                                                           \
    EXTRACTGRADIENTVARS                                                                                                \
    const Eigen::MatrixXd Iw = Eigen::MatrixXd::Identity(WV.rows(), WV.rows()); /*identity matrix of the water phase*/ \
    const Eigen::MatrixXd Ig = Eigen::MatrixXd::Identity(RH.rows(), RH.rows()); /*identity matrix of the gas phase*/   \
    /*derivatives with the naming pattern VAR1_dVAR2*/                                                                 \
    const double dcw_dwv = mDCw->d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double dcw_drh = mDCw->d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double dcg_dwv = mDCg->d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double dcg_drh = mDCg->d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double mec_dwv = mMeC->d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                          \
    const double mec_drh = mMeC->d_RelativeHumidity(wv, rh, wv_dt, rh_dt);                                             \
    const double wveq_dwv = mWVEq->d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);                                        \
    const double wveq_drh = mWVEq->d_RelativeHumidity(wv, rh, wv_dt, rh_dt);

#define EXTRACTDAMPINGVARS                                                                                             \
    EXTRACTSHAREDVARS                                                                                                  \
    const Eigen::MatrixXd GWV = Bw * WV; /*water volume graction gradient*/                                            \
    const Eigen::MatrixXd GRH = Bg * RH; /*relative humidity gradient*/                                                \
    /*derivatives with the naming pattern VAR1_dVAR2*/                                                                 \
    const double dcw_dwv_dt = mDCw->d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double dcw_drh_dt = mDCw->d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double dcg_dwv_dt = mDCg->d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double dcg_drh_dt = mDCg->d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double mec_dwv_dt = mMeC->d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                    \
    const double mec_drh_dt = mMeC->d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);                                       \
    const double wveq_dwv_dt = mWVEq->d_WaterVolumeFraction_dt(wv, rh, wv_dt, rh_dt);                                  \
    const double wveq_drh_dt = mWVEq->d_RelativeHumidity_dt(wv, rh, wv_dt, rh_dt);

namespace NuTo
{

Integrands::MoistureTransport::MoistureTransport(DofType dofTypeWV, DofType dofTypeRH,
                                                 const MTCoefficientInterface& diffCoeffWV,
                                                 const MTCoefficientInterface& diffCoeffRH,
                                                 const MTCoefficientInterface& massExchCoeff,
                                                 const MTCoefficientInterface& wvEquilibriumCoeff, double rho_w,
                                                 double rho_g_sat, double PV)
    : mDofTypeWV(dofTypeWV)
    , mDofTypeRH(dofTypeRH)
    , mRho_w(rho_w)
    , mRho_g_sat(rho_g_sat)
    , mPV(PV)
    , mDCw(diffCoeffWV.Clone())
    , mDCg(wvEquilibriumCoeff.Clone())
    , mMeC(massExchCoeff.Clone())
    , mWVEq(wvEquilibriumCoeff.Clone())
{
}

DofVector<double> Integrands::MoistureTransport::Gradient(const CellIpData& cellIpData, double deltaT)
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

DofMatrix<double> Integrands::MoistureTransport::Stiffness(const CellIpData& cellIpData, double deltaT)
{
    DofMatrix<double> stiffness;

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

DofMatrix<double> Integrands::MoistureTransport::Damping(const CellIpData& cellIpData, double deltaT)
{
    DofMatrix<double> damping;

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

DofMatrix<double> Integrands::MoistureTransport::Mass(const CellIpData& cellIpData, double deltaT)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

inline void Integrands::MoistureTransport::CheckValuesValid(double wv)
{
    if (std::abs(mPV - wv) < 10e-9)
        throw Exception(
                __PRETTY_FUNCTION__,
                "Poresvolume is completly occupied by water. Gas transport equation will produce undefined behaviour!");
}
}

#undef EXTRACTDAMPINGVARS
#undef EXTRACTSTIFFNESSVARS
#undef EXTRACTGRADIENTVARS
#undef EXTRACTSHAREDVARS
