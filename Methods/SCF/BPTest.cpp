#include "Methods/SCF/BPTest.hpp"
#include "Methods/SCF/SCF_Common.hpp"

using namespace pulsar::datastore;
using namespace pulsar::modulebase;
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

void BPTest::Initialize_(const System & sys, const std::string & bstag)
{
    const BasisSet bs = sys.GetBasisSet(bstag);

    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    mod_nuc_rep->Calculate(0, sys, &nucrep_, 1);

    ////////////////////////////
    // One-electron hamiltonian
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->SetBases(sys, bstag, bstag);
    Hcore_ = FillOneElectronMatrix(mod_ao_core, bs);


    initialized_ = true;
}


double BPTest::CalculateEnergy_(const IrrepSpinMatrixD & Dmat,
                                const IrrepSpinMatrixD & Fmat)
{
    // calculate the energy
    double energy = 0.0;
    double oneelectron = 0.0;
    double twoelectron = 0.0;

    for(auto ir : Dmat.GetIrreps())
    for(auto s : Dmat.GetSpins(ir))
    {
        const SimpleMatrixD & d = Dmat.Get(ir, s);
        const SimpleMatrixD & f = Fmat.Get(ir, s);

        for(size_t i = 0; i < d.NRows(); i++)
        for(size_t j = 0; j < d.NCols(); j++)
        {
            oneelectron += d(i,j) * Hcore_(i,j);
            twoelectron += 0.5 * d(i,j) * f(i,j);
        }
    }

    twoelectron -= 0.5*oneelectron;
    energy = oneelectron + twoelectron;

    out.Output("            One electron: %16.8e\n", oneelectron);
    out.Output("            Two electron: %16.8e\n", twoelectron);
    out.Output("        Total Electronic: %16.8e\n", energy);
    out.Output("       Nuclear Repulsion: %16.8e\n", nucrep_);

    energy += nucrep_;
    out.Output("            Total energy: %16.8e\n", energy);

    return energy;
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

    if(!initialized_)
        Initialize_(*wfn.system, bstag);

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);
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
    double damp = Options().Get<double>("DAMPING_FACTOR");


    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Load and set up the iterator  and fockbuild modules
    auto mod_iter = CreateChildFromOption<SCFIterator>("KEY_SCF_ITERATOR");
    auto mod_fock = CreateChildFromOption<FockBuilder>("KEY_FOCK_BUILDER");


    // Storing the results of the previous iterations
    // (right now, this is the initial guess)
    Wavefunction lastwfn = initial_wfn;
    double last_energy = initial_energy;  // right now, energy = initial guess energy
    double current_energy = last_energy;
    double energy_diff = 0.0;
    double dens_diff = 0.0;

    // The last density
    IrrepSpinMatrixD lastdens = FormDensity(*lastwfn.cmat, *lastwfn.occupations);
    IrrepSpinMatrixD lastfmat;



    // Start the SCF procedure
    size_t iter = 0;
    do
    {
        iter++; 

        // The Fock matrix
        IrrepSpinMatrixD Fmat = mod_fock->Build(lastwfn);

        // apply damping if we are past the first iteration
        if(iter > 1)
        {
            for(auto ir : Fmat.GetIrreps())
            for(auto s : Fmat.GetSpins(ir))
            {
                auto & m = Fmat.Get(ir, s);
                const auto & lastm = lastfmat.Get(ir, s);

                for(size_t i = 0; i < m.NRows(); i++)
                for(size_t j = 0; j < m.NCols(); j++)
                    m(i,j) = damp*lastm(i,j) + (1.0-damp)*m(i,j);

            }

        } 

        // Iterate, making a new wavefunction
        Wavefunction newwfn = mod_iter->Next(lastwfn, Fmat);

        // Form the new density and calculate the energy
        if(!newwfn.opdm)
            throw GeneralException("Returned wfn doesn't have opdm");

        const IrrepSpinMatrixD dens = *newwfn.opdm;
        current_energy = CalculateEnergy_(dens, Fmat);

        // store the energy for next time
        energy_diff = current_energy - last_energy;
        last_energy = current_energy;

        // Store this wfn as the last one
        lastwfn = newwfn;

        // Note - we've already set lastwfn equal to the new iteration
        dens_diff = CalculateRMSDens(dens, lastdens);

        // store the new density
        lastdens = std::move(dens); 
        lastfmat = std::move(Fmat);


        out.Output("%5?  %16.8e  %16.8e  %16.8e\n",
                    iter, current_energy, energy_diff, dens_diff);

    } while(fabs(energy_diff) > etol &&
            dens_diff > dtol &&
            iter < maxniter);

    //! \todo form C if only opdm is set in final wfn?

    // What are we returning 
    return {std::move(lastwfn), {current_energy}};
}
    

} // close namespace pulsarmethods
