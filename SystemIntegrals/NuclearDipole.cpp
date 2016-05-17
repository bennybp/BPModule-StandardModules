#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include <pulsar/math/PointManipulation.hpp>
#include "NuclearDipole.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::math;


uint64_t NuclearDipole::Calculate_(uint64_t deriv, const System & sys,
                                   double * outbuffer, size_t bufsize)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Nuclear Dipole with deriv != 0");

    if(bufsize < 3)
        throw GeneralException("Not enough space in output buffer");


    Point dip = WeightedPointsCenter<Point>(sys, [](const Atom & atom){ return atom.GetZ(); });
    outbuffer[0] = dip[0];
    outbuffer[1] = dip[1];
    outbuffer[2] = dip[2];

    return 1;
}
