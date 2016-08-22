#include <pulsar/modulebase/All.hpp>
#include <pulsar/util/Format.hpp> // for format_string

#include "Methods/SCF/DIIS.hpp"
#include "Methods/SCF/SCFCommon.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace pulsar::datastore;
using namespace pulsar::modulebase;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::util;
using namespace pulsar::exception;
using namespace bphash;


namespace pulsarmethods {

void DIIS::initialize_(const Wavefunction & wfn)
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


    /////////////////////// 
    // Overlap
    /////////////////////// 
    const std::string ao_overlap_key = options().get<std::string>("KEY_AO_OVERLAP");
    auto overlapimpl = mod_ao_cache->calculate(ao_overlap_key, 0, wfn, bs, bs);
    S_ = convert_to_eigen(overlapimpl.at(0));  // .at(0) = first (and only) component


    ////////////////////////////
    // One-electron hamiltonian
    ////////////////////////////
    const std::string ao_build_key = options().get<std::string>("KEY_AO_COREBUILD");
    auto Hcoreimpl = mod_ao_cache->calculate(ao_build_key, 0, wfn, bs, bs);
    Hcore_ = convert_to_eigen(Hcoreimpl.at(0));  // .at(0) = first (and only) component

    bs.print(out);
}


DIIS::DerivReturnType DIIS::deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw GeneralException("System is not set!");

    // have we already calculated this (and the result is in the cache)?
    auto hash = make_hash(HashType::Hash128, wfn);
    std::string hashstr = format_string("deriv_%?_wfn:%?", order, hash_to_string(hash));
    out.debug("Checking for key %? in the cache\n", hashstr);

    if(cache().count(hashstr))
    {
        out.debug("Found. Returning that\n");
        return *(cache().get<DerivReturnType>(hashstr));
    }

    out.debug("Not found. I have to calculate it :(\n");

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
            throw GeneralException("Missing initial guess module when I don't have a C-matrix");

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
            throw GeneralException("Missing Occupations in given wavefunction");
        if(!wfn.epsilon)
            throw GeneralException("Missing Epsilon in given wavefunction");

        initial_wfn = wfn;
    }


    ///////////////////////////////
    // Obtain the options for SCF
    ///////////////////////////////
    double etol = options().get<double>("EGY_TOLERANCE");
    size_t maxniter = options().get<size_t>("MAX_ITER");
    double dtol = options().get<double>("DENS_TOLERANCE");
    //double damp = options().get<double>("DAMPING_FACTOR");


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

    // storage for DIIS
    std::list<BlockedEigenMatrix> errqueue;
    std::list<BlockedEigenMatrix> Fqueue;

    // The last density
    IrrepSpinMatrixD lastdens = FormDensity(*lastwfn.cmat, *lastwfn.occupations);
    IrrepSpinMatrixD lastfmat;

    // for convenience
    const MatrixXd & S = *S_;

    // Start the SCF procedure
    size_t iter = 0;
    do
    {
        iter++; 

        // The Fock matrix
        IrrepSpinMatrixD Fmat = mod_fock->calculate(lastwfn);

        // pop the last one if we are beyond our limit
        if(errqueue.size() > 5)
        {
            errqueue.pop_front();
            Fqueue.pop_front();
        }

        // The error and F matrices to add
        BlockedEigenMatrix it_err;
        BlockedEigenMatrix it_f;

        // calculate the error matrix
        for(auto ir : Fmat.get_irreps())
        for(auto s : Fmat.get_spins(ir))
        {
            const MatrixXd & f = *(convert_to_eigen(Fmat.get(ir, s)));
            const MatrixXd & d = *(convert_to_eigen(lastwfn.opdm->get(ir, s)));

            MatrixXd e = f*d*S - S*d*f;
            it_err.set(ir, s, std::move(e));
            it_f.set(ir, s, f);
        }

        // add to the queues
        errqueue.push_back(std::move(it_err));
        Fqueue.push_back(std::move(it_f));

        // form the B matrix and extrapolate
        // for the new F matrix
        size_t nvec = errqueue.size();
        if(nvec > 2)
        {
            const auto & first = errqueue.front();

            for(auto ir : first.get_irreps())
            for(auto s : first.get_spins(ir))
            {
                Eigen::MatrixXd b(nvec+1, nvec+1);

                size_t i = 0;
                for(const auto & it1 : errqueue)
                {
                    size_t j = 0;
                    for(const auto & it2 : errqueue)
                    {
                        b(i,j) = (it1.get(ir, s).cwiseProduct(it2.get(ir,s))).sum();
                        j++;
                    }

                    i++;
                }

                Eigen::VectorXd r(nvec+1);

                // add the last row and column, and fill in the RHS vector
                for(i = 0; i <= nvec; i++)
                {
                    b(nvec, i) = b(i, nvec) = -1.0;
                    r(i) = 0;
                }
                b(nvec, nvec) = 0.0;
                r(nvec) = -1.0;

                // solve the system of linear equations
                Eigen::VectorXd c = b.inverse()*r;

                // form the new F matrix, and place it where we
                // got the original
                const size_t nrow = Fmat.get(ir, s)->size(0);
                const size_t ncol = Fmat.get(ir, s)->size(1);
                MatrixXd f = MatrixXd::Zero(nrow, ncol);
                i = 0;

                for(const auto & it : Fqueue)
                {
                    const MatrixXd & itf = it.get(ir, s);
                    f += c(i) * itf;
                    i++;
                }

                Fmat.set(ir, s, std::make_shared<EigenMatrixImpl>(std::move(f)));
            }
        }
        //else if(nvec == 1)
        //{
        // Do damping
        //}

        // Fmat should now have the extrapolated fock matrices



        // Iterate, making a new wavefunction
        Wavefunction newwfn = mod_iter->next(lastwfn, Fmat);


        /*
        out << "Occupations:\n";
        const auto & occ = *newwfn.occupations;
        for(auto ir : occ.get_irreps())
        {
            for(auto s : occ.get_spins(ir))
            {
                const auto & o = occ.get(ir, s);
                out << "   Spin " << s << ":";
                for(size_t q = 0; q < o.Size(); q++)
                    out << " " << o(q);
                out << "\n";
            }
        }
        */


        // Form the new density and calculate the energy
        if(!newwfn.opdm)
            throw GeneralException("Returned wfn doesn't have opdm");

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

    // cache the result
    DerivReturnType ret{std::move(lastwfn), {current_energy}};
    cache().set(hashstr, ret, CacheData::CheckpointGlobal);

    return ret;
}
    

} // close namespace pulsarmethods
