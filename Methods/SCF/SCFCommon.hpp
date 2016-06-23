#ifndef PULSAR_GUARD_SCF__SCFCOMMON_HPP_
#define PULSAR_GUARD_SCF__SCFCOMMON_HPP_

#include "Common/EigenCommon.hpp"

#include <pulsar/output/OutputStream.hpp>
#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/modulebase/TwoElectronIntegral.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>

#define INDEX2(i,j) (  (j > i) ? (j*(j+1))/2 + i : (i*(i+1))/2 + j )
#define INDEX4(i,j,k,l)  ( INDEX2(k,l) > INDEX2(i,j) ?  (INDEX2(k,l)*(INDEX2(k,l)+1))/2 + INDEX2(i,j) : (INDEX2(i,j)*(INDEX2(i,j)+1))/2 + INDEX2(k,l) )

namespace pulsarmethods {

std::vector<double>
FillTwoElectronVector(pulsar::modulemanager::ModulePtr<pulsar::modulebase::TwoElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

pulsar::math::IrrepSpinVectorD FindOccupations(size_t nelec);

pulsar::math::IrrepSpinMatrixD
FormDensity(const pulsar::math::IrrepSpinMatrixD & Cmat, const pulsar::math::IrrepSpinVectorD & occ);

double CalculateRMSDens(const pulsar::math::IrrepSpinMatrixD & m1,
                        const pulsar::math::IrrepSpinMatrixD & m2);

double Calculateenergy(const Eigen::MatrixXd & Hcore, double nucrep,
                       const pulsar::math::IrrepSpinMatrixD & Dmat,
                       const pulsar::math::IrrepSpinMatrixD & Fmat,
                       pulsar::output::OutputStream & out);

Eigen::MatrixXd FormS12(const Eigen::MatrixXd & S);

} // close namespace pulsarmethods

#endif

