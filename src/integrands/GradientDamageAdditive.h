#pragma once

#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/constitutive/LinearElasticDamage.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "nuto/mechanics/constitutive/ModifiedMisesStrainNorm.h"

namespace NuTo
{
namespace NonlocalInteraction
{
struct Constant;
}

namespace Integrands
{
//! Implicit gradient enhanced damage model
//! Peerlings RHJ et al.
//! https://dx.doi.org/10.1002/(SICI)1097-0207(19961015)39:19<3391::AID-NME7>3.0.CO;2-D
//! @tparam TDim global dimension
//! @tparam TInteraction interaction law that provides a damage dependent nonlocal parameter factor
//! @tparam TDamageLaw damage law that provides .Damage(double) and .Derivative(double)
template <int TDim, class TAdditiveIOLaw, typename TInteraction = NonlocalInteraction::Constant,
          typename TDamageLaw = Constitutive::DamageLawExponential>
class GradientDamage
{
public:
    //! ctor
    //! @param disp dof type associated with displacements
    //! @param eeq scalar dof type associated with nonlocal equivalent strains
    //! @param c nonlocal parameter unit length squared
    //! @param linearElasticDamage damage law that avoids adding damage to the negative volumetric part
    //! @param additiveIOLaw: Additional law that affect the output stresses and the input strains
    //! @param damageLaw damage law that provides .Damage(double) and .Derivative(double)
    //! @param strainNorm modified mises strain norm
    //! @param interaction damage dependent nonlocal parameter factor
    GradientDamage(DofType disp, ScalarDofType eeq, double c, Laws::LinearElasticDamage<TDim> linearElasticDamage,
                   TAdditiveIOLaw additiveIOLaw, TDamageLaw damageLaw,
                   Constitutive::ModifiedMisesStrainNorm<TDim> strainNorm, TInteraction interaction = TInteraction())
        : mDisp(disp)
        , mEeq(eeq)
        , mC(c)
        , mLinearElasticDamage(linearElasticDamage)
        , mAdditiveIOLaw{additiveIOLaw}
        , mDamageLaw(damageLaw)
        , mNorm(strainNorm)
        , mInteraction(interaction)
    {
    }

    //! ctor
    //! @param disp dof type associated with displacements
    //! @param eeq scalar dof type associated with nonlocal equivalent strains
    //! @param m softening material parameter
    //! @param additiveIOLaw: Additional law that affect the output stresses and the input strains
    //! @param damageApplication option to skip damaging the hydrostatic compressive part
    //! @param interaction damage dependent nonlocal parameter factor
    GradientDamage(DofType disp, ScalarDofType eeq, Material::Softening m, TAdditiveIOLaw additiveIOLaw,
                   Laws::eDamageApplication damageApplication = Laws::eDamageApplication::FULL,
                   TInteraction interaction = TInteraction())
        : mDisp(disp)
        , mEeq(eeq)
        , mC(m.c)
        , mLinearElasticDamage(m, damageApplication)
        , mAdditiveIOLaw{additiveIOLaw}
        , mDamageLaw(m)
        , mNorm(m)
        , mInteraction(interaction)
    {
    }

    virtual ~GradientDamage() = default;

    DofVector<double> Gradient(const CellIpData& data)
    {
        DofVector<double> gradient;

        double eeq = data.Value(mEeq);
        double omega = mDamageLaw.Damage(Kappa(data));

        auto eeqGradient = data.Apply(mEeq, Nabla::Gradient());
        NuTo::EngineeringStrain<TDim> strainMechanics = StrainMechanics(data);

        Eigen::MatrixXd Neeq = data.N(mEeq);
        Eigen::MatrixXd Beeq = data.B(mEeq, Nabla::Gradient());
        Eigen::MatrixXd Bdisp = data.B(mDisp, Nabla::Strain());

        double g = mInteraction.Factor(omega);

        gradient[mDisp] = Bdisp.transpose() * Stress(data);
        gradient[mEeq] =
                Neeq.transpose() * (eeq - mNorm.Value(strainMechanics)) + Beeq.transpose() * (mC * g * eeqGradient);

        return gradient;
    }

    DofMatrix<double> Hessian0(const CellIpData& data)
    {
        DofMatrix<double> hessian0;

        double kappa = Kappa(data);
        double omega = mDamageLaw.Damage(kappa);
        double dKappa_dEeq = DkappaDeeq(data);
        double dOmega_dKappa = mDamageLaw.Derivative(kappa);

        NuTo::EngineeringStrain<TDim> strainMechanics = StrainMechanics(data);

        Eigen::MatrixXd Neeq = data.N(mEeq);
        Eigen::MatrixXd Beeq = data.B(mEeq, Nabla::Gradient());
        Eigen::MatrixXd Bdisp = data.B(mDisp, Nabla::Strain());

        double g = mInteraction.Factor(omega);
        double dgdw = mInteraction.Derivative(omega);
        auto eeqGradient = data.Apply(mEeq, Nabla::Gradient());

        hessian0(mDisp, mDisp) =
                Bdisp.transpose() * mLinearElasticDamage.DstressDstrain(strainMechanics, omega) * Bdisp;
        hessian0(mEeq, mDisp) = -Neeq.transpose() * mNorm.Derivative(strainMechanics).transpose() * Bdisp;
        hessian0(mEeq, mEeq) = Neeq.transpose() * Neeq + mC * g * Beeq.transpose() * Beeq +
                               Beeq.transpose() * mC * eeqGradient * dgdw * dOmega_dKappa * dKappa_dEeq * Neeq;
        hessian0(mDisp, mEeq) =
                Bdisp.transpose() *
                (mLinearElasticDamage.DstressDomega(strainMechanics, omega) * dOmega_dKappa * dKappa_dEeq) * Neeq;

        return hessian0;
    }


    Eigen::VectorXd Strain(const CellIpData& data) const
    {
        return data.Apply(mDisp, Nabla::Strain());
    }


    Eigen::VectorXd StrainMechanics(const CellIpData& data) const
    {
        return Strain(data) - StrainAdditive(data);
    }

    Eigen::VectorXd StrainAdditive(const CellIpData& data) const
    {
        return mAdditiveIOLaw.Strain(data, 0);
    }


    Eigen::VectorXd Stress(const CellIpData& data) const
    {
        return StressMechanics(data) + StressAdditive(data);
    }

    Eigen::VectorXd StressMechanics(const CellIpData& data) const
    {
        double omega = mDamageLaw.Damage(Kappa(data));
        return mLinearElasticDamage.Stress(StrainMechanics(data), omega);
    }

    Eigen::VectorXd StressAdditive(const CellIpData& data) const
    {
        return mAdditiveIOLaw.Stress(data, 0.);
    }

    virtual void Update(const CellIpData& data)
    {
        mKappas(data.Ids().cellId, data.Ids().ipId) = Kappa(data);
    }

    virtual double Kappa(const CellIpData& data) const
    {
        return std::max(mKappas(data.Ids().cellId, data.Ids().ipId), data.Value(mEeq));
    }

    virtual double DkappaDeeq(const CellIpData& data) const
    {
        return data.Value(mEeq) >= mKappas(data.Ids().cellId, data.Ids().ipId) ? 1 : 0;
    }

    Eigen::MatrixXd mKappas;

    DofType mDisp;
    ScalarDofType mEeq;
    double mC;
    Laws::LinearElasticDamage<TDim> mLinearElasticDamage;
    TAdditiveIOLaw mAdditiveIOLaw;
    TDamageLaw mDamageLaw;
    Constitutive::ModifiedMisesStrainNorm<TDim> mNorm;
    TInteraction mInteraction;
};
} /* Integrand */

namespace NonlocalInteraction
{

//! Results in the model used by Peerlings et al.
struct Constant
{
    double Factor(double) const
    {
        return 1;
    }
    double Derivative(double) const
    {
        return 0;
    }
};

//! Results in the model used by Poh & Sun 2017, IJNME and limits the nonlocal parameter
struct Decreasing
{
    Decreasing(double R = 0.005, double eta = 5)
        : mR(R)
        , mEta(eta)
    {
    }

    double Factor(double omega) const
    {
        return ((1. - mR) * std::exp(-mEta * omega) + mR - std::exp(-mEta)) / (1. - std::exp(-mEta));
    }
    double Derivative(double omega) const
    {
        return ((1. - mR) * std::exp(-mEta * omega)) / (1. - std::exp(-mEta)) * -mEta;
    }
    double mR;
    double mEta;
};
} /* NonlocalInteraction */
} /* NuTo */
