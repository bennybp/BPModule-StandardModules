#ifndef BPTEST_COMMON_HPP
#define BPTEST_COMMON_HPP

#include <pulsar/modulebase/All.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/BlockByIrrepSpin.hpp>

#include <eigen3/Eigen/Dense>

#include <vector>


namespace pulsarmethods {

#define INDEX2(i,j) (  (j > i) ? (j*(j+1))/2 + i : (i*(i+1))/2 + j )
#define INDEX4(i,j,k,l)  ( INDEX2(k,l) > INDEX2(i,j) ?  (INDEX2(k,l)*(INDEX2(k,l)+1))/2 + INDEX2(i,j) : (INDEX2(i,j)*(INDEX2(i,j)+1))/2 + INDEX2(k,l) )

Eigen::MatrixXd
FillOneElectronMatrix(pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

std::vector<double>
FillTwoElectronVector(pulsar::modulemanager::ModulePtr<pulsar::modulebase::TwoElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

pulsar::math::BlockByIrrepSpin<Eigen::MatrixXd>
FormDensity(const pulsar::math::BlockByIrrepSpin<Eigen::MatrixXd> & Cmat,
            const pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> & occ);


pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> FindOccupations(size_t nelec);

pulsar::math::SimpleMatrixD EigenToSimpleMatrix(const Eigen::MatrixXd & m);

pulsar::math::SimpleVectorD EigenToSimpleVector(const Eigen::VectorXd & v);

Eigen::MatrixXd SimpleMatrixToEigen(const pulsar::math::SimpleMatrixD & m); 

Eigen::VectorXd SimpleVectorToEigen(const pulsar::math::SimpleVectorD & v); 

} // close namespace pulsarmethods

#endif

