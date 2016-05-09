#ifndef BPTEST_COMMON_HPP
#define BPTEST_COMMON_HPP

#include <pulsar/modulebase/All.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/BlockByIrrepSpin.hpp>

#include <eigen3/Eigen/Dense>

#include <vector>


namespace pulsarmethods {

Eigen::MatrixXd
FillOneElectronMatrix(pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);



pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> FindOccupations(size_t nelec);

pulsar::math::SimpleMatrixD EigenToSimpleMatrix(const Eigen::MatrixXd & m);

pulsar::math::SimpleVectorD EigenToSimpleVector(const Eigen::VectorXd & v);

Eigen::MatrixXd SimpleMatrixToEigen(const pulsar::math::SimpleMatrixD & m); 

Eigen::VectorXd SimpleVectorToEigen(const pulsar::math::SimpleVectorD & v); 

} // close namespace pulsarmethods

#endif

