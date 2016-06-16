#include "Methods/SCF/DIIS.hpp"
#include "Methods/SCF/SCF_Common.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace pulsar::datastore;
using namespace pulsar::modulebase;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::exception;


namespace pulsarmethods {

void DIIS::Initialize_(const Wavefunction & wfn)
{
    // get the basis set
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    const BasisSet bs = wfn.system->GetBasisSet(bstag);

    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    mod_nuc_rep->Initialize(0, *wfn.system);
    mod_nuc_rep->Calculate(&nucrep_, 1);


    ////////////////////////////
    // Overlap
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->Initialize(0, wfn, bs, bs);
    S_ = FillOneElectronMatrix(mod_ao_overlap, bs);


    ////////////////////////////
    // One-electron hamiltonian
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->Initialize(0, wfn, bs, bs);
    Hcore_ = FillOneElectronMatrix(mod_ao_core, bs);

    bs.Print(out);
}


double DIIS::CalculateEnergy_(const IrrepSpinMatrixD & Dmat,
                              const IrrepSpinMatrixD & Fmat)
{
    // calculate the energy
    double energy = 0.0;
    double oneelectron = 0.0;
    double twoelectron = 0.0;

    for(auto ir : Dmat.GetIrreps())
    for(auto s : Dmat.GetSpins(ir))
    {
        const MatrixXd & d = *(convert_to_eigen(Dmat.Get(ir, s)));
        const MatrixXd & f = *(convert_to_eigen(Fmat.Get(ir, s)));

        for(long i = 0; i < d.rows(); i++)
        for(long j = 0; j < d.cols(); j++)
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



DIIS::DerivReturnType DIIS::Deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw GeneralException("System is not set!");

    Initialize_(wfn); // will only use the system from the wfn

    std::string bstag = Options().Get<std::string>("BASIS_SET");

    const BasisSet bs = wfn.system->GetBasisSet(bstag);
  
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
    //double damp = Options().Get<double>("DAMPING_FACTOR");


    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Load and set up the iterator  and fockbuild modules
    auto mod_iter = CreateChildFromOption<SCFIterator>("KEY_SCF_ITERATOR");
    auto mod_fock = CreateChildFromOption<FockBuilder>("KEY_FOCK_BUILDER");
    mod_fock->Initialize(order, wfn, bs); 


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

    // Start the SCF procedure
    size_t iter = 0;
    do
    {
        iter++; 

        // The Fock matrix
        IrrepSpinMatrixD Fmat = mod_fock->Calculate(lastwfn);

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
        for(auto ir : Fmat.GetIrreps())
        for(auto s : Fmat.GetSpins(ir))
        {
            const MatrixXd & f = *(convert_to_eigen(Fmat.Get(ir, s)));
            const MatrixXd & d = *(convert_to_eigen(lastwfn.opdm->Get(ir, s)));

            MatrixXd e = f*d*S_ - S_*d*f;
            it_err.Take(ir, s, std::move(e));
            it_f.Set(ir, s, f);
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

            for(auto ir : first.GetIrreps())
            for(auto s : first.GetSpins(ir))
            {
                Eigen::MatrixXd b(nvec+1, nvec+1);

                size_t i = 0;
                for(const auto & it1 : errqueue)
                {
                    size_t j = 0;
                    for(const auto & it2 : errqueue)
                    {
                        b(i,j) = (it1.Get(ir, s).cwiseProduct(it2.Get(ir,s))).sum();
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
                const size_t nrow = Fmat.Get(ir, s)->size(0);
                const size_t ncol = Fmat.Get(ir, s)->size(1);
                MatrixXd f = MatrixXd::Zero(nrow, ncol);
                i = 0;

                for(const auto & it : Fqueue)
                {
                    const MatrixXd & itf = it.Get(ir, s);
                    f += c(i) * itf;
                    i++;
                }

                Fmat.Take(ir, s, std::make_shared<EigenMatrixImpl>(std::move(f)));
            }
        }
        //else if(nvec == 1)
        //{
        // Do damping
        //}

        // Fmat should now have the extrapolated fock matrices



        // Iterate, making a new wavefunction
        Wavefunction newwfn = mod_iter->Next(lastwfn, Fmat);


        /*
        out << "Occupations:\n";
        const auto & occ = *newwfn.occupations;
        for(auto ir : occ.GetIrreps())
        {
            for(auto s : occ.GetSpins(ir))
            {
                const auto & o = occ.Get(ir, s);
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

    } while(fabs(energy_diff) > etol ||
            dens_diff > dtol &&
            iter < maxniter);

    //! \todo form C if only opdm is set in final wfn?

    // What are we returning 
    return {std::move(lastwfn), {current_energy}};
}
    

} // close namespace pulsarmethods
