#ifndef PULSAR_GUARD_ONEELECTRONINTEGRALS__OSOVERLAP_HPP_
#define PULSAR_GUARD_ONEELECTRONINTEGRALS__OSOVERLAP_HPP_

void OSOverlap(const double alpha1, const double xyz1[3],
               const double alpha2, const double xyz2[3],
               int nam1, int nam2,
               double ** outbuffer);

#endif
