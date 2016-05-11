#include "Methods/SCF/BPTest.hpp"
#include "Methods/SCF/SCF_Common.hpp"

#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/math/Cast.hpp>

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::datastore;


namespace pulsarmethods{



BPTest::DerivReturnType BPTest::Deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);
    size_t nao = bs.NFunctions();
    size_t nshell = bs.NShell();
    size_t maxnfunc = bs.MaxNFunctions();
    size_t maxnfunc2 = maxnfunc * maxnfunc;
    bs.Print(out);
  
    //////////////////////////////////////////////////////
    // Storage of eigen matrices, etc, by irrep and spin 
    // These represent the "current" cmatrix, etc
    //////////////////////////////////////////////////////
    Wavefunction initial_wfn;
    double initial_energy = 0;
 
    //////////////////////////
    // Initial Guess
    //////////////////////////
    if(!wfn.cmat) // c-matrix hasn't been set in the passed wfn
    {
        out.Debug("Don't have C-matrices set. Will call initial guess module\n");

        if(!Options().Has("KEY_INITIAL_GUESS"))
            throw GeneralException("Missing initial guess module when I don't have a C-matrix");

        // load and run the initial guess module
        auto mod_iguess = CreateChildFromOption<EnergyMethod>("KEY_INITIAL_GUESS");
        auto iguess_ret = mod_iguess->Energy(wfn);

        initial_wfn = iguess_ret.first;
        initial_energy = iguess_ret.second;
    }
    else
    {
        out.Debug("Using initial wavefunction as a starting point");
 
        // c-matrices have been set. Make sure we have occupations, etc, as well
        if(!wfn.occupations)
            throw GeneralException("Missing Occupations");
        if(!wfn.epsilon)
            throw GeneralException("Missing Epsilon");

        initial_wfn = wfn;
    }


    ///////////////////////////////
    // Obtain the options for SCF
    ///////////////////////////////
    double etol = Options().Get<double>("E_TOLERANCE");
    size_t maxniter = Options().Get<double>("MAX_ITER");


    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Storing the results of the previous iterations
    // (right now, this is the initial guess)
    Wavefunction lastwfn = initial_wfn;
    double last_energy = initial_energy;  // right now, energy = initial guess energy
    double current_energy = last_energy;
    double energy_diff = 0.0;

    // Start the SCF procedure
    size_t iter = 0;
    do
    {
        iter++; 

        auto mod_iter = CreateChildFromOption<EnergyMethod>("KEY_SCF_ITERATOR");

        // Iterate, making a new wavefunction
        auto iterate_ret = mod_iter->Energy(lastwfn);

        // New energy
        last_energy = current_energy;
        current_energy = iterate_ret.second;
        lastwfn = iterate_ret.first;

        out.Output("Iteration %?\n", iter);
        out.Output("       Total energy: %16.8e\n", current_energy);

        energy_diff = current_energy - last_energy;

        out.Output(" Difference from last step: %16.8e\n", energy_diff);

    } while(fabs(energy_diff) > etol && iter < maxniter);

    // What are we returning 
    return {std::move(lastwfn), {current_energy}};
}
    

}//End namespace
