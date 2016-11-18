#include <pulsar/modulebase/All.hpp>
#include "Methods/SCF/Damping.hpp"
#include "Methods/SCF/SCFCommon.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace pulsar;


namespace pulsarmethods {

void Damping::initialize_(const Wavefunction & wfn)
{
    // get the basis set
    std::string bstag = options().get<std::string>("BASIS_SET");

    const BasisSet bs = wfn.system->get_basis_set(bstag);

    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = create_child_from_option<SystemIntegral>("KEY_NUC_REPULSION");
    mod_nuc_rep->initialize(0, *wfn.system);
    mod_nuc_rep->calculate(&nucrep_, 1);

    /////////////////////////////////////
    // The one-electron integral cacher
    /////////////////////////////////////
    auto mod_ao_cache = create_child_from_option<OneElectronMatrix>("KEY_ONEEL_MAT");

    ////////////////////////////
    // One-electron hamiltonian
    ////////////////////////////
    const std::string ao_build_key = options().get<std::string>("KEY_AO_COREBUILD");
    auto Hcoreimpl = mod_ao_cache->calculate(ao_build_key, 0, wfn, bs, bs);
    Hcore_ = convert_to_eigen(*Hcoreimpl.at(0));  // .at(0) = first (and only) component

    bs.print(out);
}


DerivReturnType Damping::deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw PulsarException("System is not set!");

    initialize_(wfn); // will only use the system from the wfn

    std::string bstag = options().get<std::string>("BASIS_SET");
    const BasisSet bs = wfn.system->get_basis_set(bstag);
  
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
        out.debug("Don't have C-matrices set. Will call initial guess module\n");

        if(!options().has("KEY_INITIAL_GUESS"))
            throw PulsarException("Missing initial guess module when I don't have a C-matrix");

        // load and run the initial guess module
        auto mod_iguess = create_child_from_option<EnergyMethod>("KEY_INITIAL_GUESS");
        auto iguess_ret = mod_iguess->energy(wfn);

        initial_wfn = iguess_ret.first;
        initial_energy = iguess_ret.second;
    }
    else
    {
        out.debug("Using given wavefunction as a starting point");
 
        // c-matrices have been set. Make sure we have occupations, etc, as well
        if(!wfn.occupations)
            throw PulsarException("Missing Occupations in given wavefunction");
        if(!wfn.epsilon)
            throw PulsarException("Missing Epsilon in given wavefunction");

        initial_wfn = wfn;
    }


    ///////////////////////////////
    // Obtain the options for SCF
    ///////////////////////////////
    double etol = options().get<double>("EGY_TOLERANCE");
    size_t maxniter = options().get<size_t>("MAX_ITER");
    double dtol = options().get<double>("DENS_TOLERANCE");
    double damp = options().get<double>("DAMPING_FACTOR");


    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Load and set up the iterator  and fockbuild modules
    auto mod_iter = create_child_from_option<SCFIterator>("KEY_SCF_ITERATOR");
    auto mod_fock = create_child_from_option<FockBuilder>("KEY_FOCK_BUILDER");
    mod_fock->initialize(order, wfn, bs); 


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
        IrrepSpinMatrixD Fmat = mod_fock->calculate(lastwfn);

        // apply damping if we are past the first iteration
        if(iter > 1)
        {
            for(auto ir : Fmat.get_irreps())
            for(auto s : Fmat.get_spins(ir))
            {
                // making a copy
                MatrixXd m = *(convert_to_eigen(*Fmat.get(ir, s)));
                const MatrixXd & lastm = *(convert_to_eigen(*lastfmat.get(ir, s)));

                for(long i = 0; i < m.rows(); i++)
                for(long j = 0; j < m.cols(); j++)
                    m(i,j) = damp*lastm(i,j) + (1.0-damp)*m(i,j);

                Fmat.set(ir, s, std::make_shared<EigenMatrixImpl>(std::move(m)));

            }
        } 

        // Iterate, making a new wavefunction
        Wavefunction newwfn = mod_iter->next(lastwfn, Fmat);

        // Form the new density and calculate the energy
        if(!newwfn.opdm)
            throw PulsarException("Returned wfn doesn't have opdm");

        const IrrepSpinMatrixD dens = *newwfn.opdm;
        current_energy = Calculateenergy(*Hcore_, nucrep_, dens, Fmat, out);

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


        out.output("%5?  %16.8e  %16.8e  %16.8e\n",
                    iter, current_energy, energy_diff, dens_diff);

    } while(fabs(energy_diff) > etol ||
            dens_diff > dtol &&
            iter < maxniter);

    //! \todo form C if only opdm is set in final wfn?

    // What are we returning 
    return {std::move(lastwfn), {current_energy}};
}
    

} // close namespace pulsarmethods
