#ifndef PULSAR_GUARD_SCF__SCF_COMMON_HPP_
#define PULSAR_GUARD_SCF__SCF_COMMON_HPP_

#include <pulsar/modulebase/All.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/BlockByIrrepSpin.hpp>
#include <pulsar/datastore/Wavefunction.hpp>

#include <eigen3/Eigen/Dense>

#include <vector>


namespace pulsarmethods {

typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MappedMatrix;
typedef Eigen::Map<Eigen::VectorXd> MappedVector;
typedef Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MappedConstMatrix;
typedef Eigen::Map<const Eigen::VectorXd> MappedConstVector;

typedef pulsar::math::BlockByIrrepSpin<Eigen::MatrixXd> BlockedEigenMatrix;
typedef pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> BlockedEigenVector;


#define INDEX2(i,j) (  (j > i) ? (j*(j+1))/2 + i : (i*(i+1))/2 + j )
#define INDEX4(i,j,k,l)  ( INDEX2(k,l) > INDEX2(i,j) ?  (INDEX2(k,l)*(INDEX2(k,l)+1))/2 + INDEX2(i,j) : (INDEX2(i,j)*(INDEX2(i,j)+1))/2 + INDEX2(k,l) )

Eigen::MatrixXd
FillOneElectronMatrix(pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

std::vector<double>
FillTwoElectronVector(pulsar::modulemanager::ModulePtr<pulsar::modulebase::TwoElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

pulsar::math::IrrepSpinMatrixD
FormDensity(const pulsar::math::IrrepSpinMatrixD & Cmat, const pulsar::math::IrrepSpinVectorD & occ);


pulsar::math::SimpleMatrixD
FormDensity(const pulsar::math::SimpleMatrixD & Cmat,
            const pulsar::math::SimpleVectorD & occ);

pulsar::math::IrrepSpinVectorD FindOccupations(size_t nelec);

pulsar::math::SimpleMatrixD EigenToSimpleMatrix(const Eigen::MatrixXd & m);

pulsar::math::SimpleVectorD EigenToSimpleVector(const Eigen::VectorXd & v);

Eigen::MatrixXd SimpleMatrixToEigen(const pulsar::math::SimpleMatrixD & m); 

Eigen::VectorXd SimpleVectorToEigen(const pulsar::math::SimpleVectorD & v); 

MappedMatrix MapSimpleMatrix(pulsar::math::SimpleMatrixD & m);

MappedConstMatrix MapConstSimpleMatrix(const pulsar::math::SimpleMatrixD & m);

MappedVector MapSimpleVector(pulsar::math::SimpleVectorD & v);

MappedConstVector MapConstSimpleVector(const pulsar::math::SimpleVectorD & v);

double CalculateRMSDens(const pulsar::math::IrrepSpinMatrixD & m1,
                        const pulsar::math::IrrepSpinMatrixD & m2);


} // close namespace pulsarmethods

#endif

