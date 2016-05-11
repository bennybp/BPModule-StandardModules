#include "Methods/SCF/BPTest.hpp"
#include "Methods/SCF/SCF_Common.hpp"

using namespace pulsar::datastore;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::exception;


namespace pulsarmethods {

static double CalculateRMSDens(const IrrepSpinMatrixD & m1, const IrrepSpinMatrixD & m2)
{
    if(!m1.SameStructure(m2))
        throw GeneralException("Density matrices have different structure");

    double rms = 0.0;

    for(Irrep ir : m1.GetIrreps())
    for(int spin : m1.GetSpins(ir))
    {
        const auto & mat1 = m1.Get(ir, spin);
        const auto & mat2 = m2.Get(ir, spin);

        for(size_t i = 0; i < mat1.NRows(); i++)
        for(size_t j = 0; j < mat1.NCols(); j++)
        {
            const double diff = mat1(i,j) - mat2(i,j);
            rms += diff*diff;
        }
    }

    return sqrt(rms);
}



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
        out.Debug("Using given wavefunction as a starting point");
 
        // c-matrices have been set. Make sure we have occupations, etc, as well
        if(!wfn.occupations)
            throw GeneralException("Missing Occupations in given wavefunction");
        if(!wfn.epsilon)
            throw GeneralException("Missing Epsilon in given wavefunction");

        initial_wfn = wfn;
    }


    ///////////////////////////////
    // Obtain the options for SCF
    ///////////////////////////////
    double etol = Options().Get<double>("E_TOLERANCE");
    size_t maxniter = Options().Get<size_t>("MAX_ITER");
    double dtol = Options().Get<double>("DENS_TOLERANCE");


    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Storing the results of the previous iterations
    // (right now, this is the initial guess)
    Wavefunction lastwfn = initial_wfn;
    double last_energy = initial_energy;  // right now, energy = initial guess energy
    double current_energy = last_energy;
    double energy_diff = 0.0;
    double dens_diff = 0.0;

    // The last density
    IrrepSpinMatrixD lastdens = FormDensity(*lastwfn.cmat, *lastwfn.occupations);

    // Load and set up the iterator module
    auto mod_iter = CreateChildFromOption<EnergyMethod>("KEY_SCF_ITERATOR");

    // Start the SCF procedure
    size_t iter = 0;
    do
    {
        iter++; 

        // Iterate, making a new wavefunction
        auto iterate_ret = mod_iter->Energy(lastwfn);

        // New energy
        last_energy = current_energy;
        current_energy = iterate_ret.second;
        lastwfn = iterate_ret.first;

        energy_diff = current_energy - last_energy;

        // Form the density and calculate the RMS difference
        // Note - we've already set lastwfn equal to the new iteration
        IrrepSpinMatrixD dens = FormDensity(*lastwfn.cmat, *lastwfn.occupations);
        dens_diff = CalculateRMSDens(dens, lastdens);

        // store the new density
        lastdens = std::move(dens); 

        out.Output("%5?  %16.8e  %16.8e  %16.8e\n",
                    iter, current_energy, energy_diff, dens_diff);

    } while(fabs(energy_diff) > etol &&
            dens_diff > dtol &&
            iter < maxniter);

    // What are we returning 
    return {std::move(lastwfn), {current_energy}};
}
    

} // close namespace pulsarmethods
