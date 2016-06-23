#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>

#include "Integrals/NuclearRepulsion.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;


void NuclearRepulsion::initialize_(unsigned int deriv, const System & sys)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Nuclear Dipole with deriv != 0");

    sys_ = &sys;
}

uint64_t NuclearRepulsion::calculate_(double * outbuffer, size_t bufsize)
{
    if(bufsize == 0)
        throw GeneralException("Not enough space in output buffer");


    double enuc = 0.0;
    for(auto it1 = sys_->begin(); it1 != sys_->end(); ++it1)
    {
        auto it2 = it1;
        std::advance(it2, 1);

        for(; it2 != sys_->end(); ++it2)
            enuc += (it1->Z * it2->Z) / it1->distance(*it2);
    }

    outbuffer[0] = enuc;
    return 1;
}
