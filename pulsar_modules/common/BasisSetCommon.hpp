#ifndef _GUARD_INTEGRAL_COMMON_HPP_
#define _GUARD_INTEGRAL_COMMON_HPP_

#include <pulsar/system/BasisSet.hpp>
#include <pulsar/datastore/CacheData.hpp>
#include <pulsar/output/OutputStream.hpp>


std::shared_ptr<const pulsar::BasisSet>
NormalizeBasis(pulsar::CacheData & cache,
               pulsar::OutputStream & out,
               const pulsar::BasisSet & bs);



#endif
