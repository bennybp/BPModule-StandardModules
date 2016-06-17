#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include <pulsar/math/PointManipulation.hpp>

#include "Integrals/NuclearDipole.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::math;


void NuclearDipole::Initialize_(unsigned int deriv, const System & sys)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Nuclear Dipole with deriv != 0");

    sys_ = &sys;
}

uint64_t NuclearDipole::Calculate_(double * outbuffer, size_t bufsize)
{
    if(bufsize < 3)
        throw GeneralException("Not enough space in output buffer");

    outbuffer[0] = 0.0;
    outbuffer[1] = 0.0;
    outbuffer[2] = 0.0;

    for(const auto & atom : *sys_)
    {
        CoordType c = atom.GetCoords();
        outbuffer[0] += atom.GetZ()*c[0];
        outbuffer[1] += atom.GetZ()*c[1];
        outbuffer[2] += atom.GetZ()*c[2];
    }

    return 3;
}
