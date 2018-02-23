#pragma once
#include "MoistureTransport.h"


template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::Integrands::MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::MoistureTransport(DofType dofTypeRH,
                                                                                      DofType dofTypeWV, double rho_w,
                                                                                      double rho_g_sat, double PV)
    : mDofTypeRH(dofTypeRH)
    , mDofTypeWV(dofTypeWV)
    , mRho_w(rho_w)
    , mRho_g_sat(rho_g_sat)
    , mPV(PV)
{
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofVector<double> NuTo::Integrands::MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Gradient(
        const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double deltaT)
{
    DofVector<double> gradient;

    Eigen::VectorXd WV = cellData.GetNodeValues(mDofTypeWV);
    Eigen::MatrixXd Nw = cellIpData.GetNMatrix(mDofTypeWV);
    Eigen::MatrixXd Bw = cellIpData.GetBMatrixGradient(mDofTypeWV);
    Eigen::VectorXd RH = cellData.GetNodeValues(mDofTypeRH);
    Eigen::MatrixXd Ng = cellIpData.GetNMatrix(mDofTypeRH);
    Eigen::MatrixXd Bg = cellIpData.GetBMatrixGradient(mDofTypeRH);

    double wv = (Nw * WV)[0];
    double wv_dt = 0.0;
    double rh = (Ng * RH)[0];
    double rh_dt = 0.0;

    double dcw = TDCw::value(wv, rh, wv_dt, rh_dt);
    double dcg = TDCg::value(wv, rh, wv_dt, rh_dt);
    double mec = TMeC::value(wv, rh, wv_dt, rh_dt);
    double wveq = TWVEq::value(wv, rh, wv_dt, rh_dt);
    double me = (wveq - wv) * mec;

    // clang-format off
    gradient[mDofTypeWV] = -dcw * Bw.transpose() * Bw * WV
                         + (me - mRho_w * wv_dt) * Nw.transpose();

    gradient[mDofTypeRH] = -dcg * Bg.transpose() * Bg * RH
                         + (mRho_g_sat * ((wv - mPV) * rh_dt + rh * wv_dt) - me) * Ng.transpose();
    // clang-format on
    return gradient;
}

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofMatrix<double> NuTo::Integrands::MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Stiffness(
        const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double deltaT)
{
    NuTo::DofMatrix<double> stiffness;

    Eigen::VectorXd WV = cellData.GetNodeValues(mDofTypeWV);
    Eigen::MatrixXd Nw = cellIpData.GetNMatrix(mDofTypeWV);
    Eigen::MatrixXd Bw = cellIpData.GetBMatrixGradient(mDofTypeWV);
    Eigen::VectorXd RH = cellData.GetNodeValues(mDofTypeRH);
    Eigen::MatrixXd Ng = cellIpData.GetNMatrix(mDofTypeRH);
    Eigen::MatrixXd Bg = cellIpData.GetBMatrixGradient(mDofTypeRH);

    double wv = (Nw * WV)[0];
    double wv_dt = 0.0;
    double rh = (Ng * RH)[0];
    double rh_dt = 0.0;

    double dcw = TDCw::value(wv, rh, wv_dt, rh_dt);
    double dcg = TDCg::value(wv, rh, wv_dt, rh_dt);
    double mec = TMeC::value(wv, rh, wv_dt, rh_dt);
    double wveq = TWVEq::value(wv, rh, wv_dt, rh_dt);

    int numDofsWV = WV.rows();
    int numDofsRH = RH.rows();
    Eigen::MatrixXd Iw = Eigen::MatrixXd::Identity(numDofsWV, numDofsWV);
    Eigen::MatrixXd Ig = Eigen::MatrixXd::Identity(numDofsRH, numDofsRH);
    double dcw_dwv = TDCw::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);
    double dcw_drh = TDCw::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);
    double dcg_dwv = TDCg::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);
    double dcg_drh = TDCg::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);
    double mec_dwv = TMeC::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);
    double mec_drh = TMeC::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);
    double wveq_dwv = TWVEq::d_WaterVolumeFraction(wv, rh, wv_dt, rh_dt);
    double wveq_drh = TWVEq::d_RelativeHumidity(wv, rh, wv_dt, rh_dt);

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

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
NuTo::DofMatrix<double> NuTo::Integrands::MoistureTransport<TDim, TDCw, TDCg, TMeC, TWVEq>::Damping(
        const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double deltaT)
{
    NuTo::DofMatrix<double> damping;

    Eigen::VectorXd WV = cellData.GetNodeValues(mDofTypeWV);
    Eigen::MatrixXd Nw = cellIpData.GetNMatrix(mDofTypeWV);
    Eigen::MatrixXd Bw = cellIpData.GetBMatrixGradient(mDofTypeWV);
    Eigen::VectorXd RH = cellData.GetNodeValues(mDofTypeRH);
    Eigen::MatrixXd Ng = cellIpData.GetNMatrix(mDofTypeRH);
    Eigen::MatrixXd Bg = cellIpData.GetBMatrixGradient(mDofTypeRH);

    damping(mDofTypeWV, mDofTypeWV) = -Bw.transpose() * Bw;
    damping(mDofTypeWV, mDofTypeRH) = -Bw.transpose() * Bg;
    damping(mDofTypeRH, mDofTypeWV) = -Bg.transpose() * Bw;
    damping(mDofTypeRH, mDofTypeRH) = -Bg.transpose() * Bg;

    return damping;
}


// template class NuTo::Integrands::MoistureTransport<1>;
// template class NuTo::Integrands::MoistureTransport<2>;
// template class NuTo::Integrands::MoistureTransport<3>;
