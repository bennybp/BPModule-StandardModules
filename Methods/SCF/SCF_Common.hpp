#ifndef BPTEST_COMMON_HPP
#define BPTEST_COMMON_HPP

#include <pulsar/modulebase/All.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>

#include <eigen3/Eigen/Dense>

#include <vector>


namespace pulsarmethods {

Eigen::MatrixXd
FillOneElectronMatrix(pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

}

#endif

