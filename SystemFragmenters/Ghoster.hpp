/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/** \file Ghoster.hpp
 *  \brief Write Me!!!!!!
 *  \author Ryan M. Richard
 *  \version 1.0
 *  \date June 6, 2016
 */

#ifndef PULSAR_GUARD_GHOSTER_HPP
#define PULSAR_GUARD_GHOSTER_HPP

#include<pulsar/modulebase/SystemFragmenter.hpp>

/** \brief This class fragments a molecule and then applies ghost atoms in a 
 *         specified manner
 * 
 *  There are really two ghosting schemes we are interested in: the Boys and
 *  Bernardi scheme where each fragment is given the rest of the system as
 *  ghosts and Valiron-Mayer Functional counterpoise where for an N-body MBE
 *  truncated at order \f$n\f$, we need to consider unions of real 
 *  \f$r\$-mers and ghost \f$g\$-mers such
 *  that \f$r+g\le n\f$.  Also of some interest is the many-ghost, many-body
 *  expansion which allows you finer control of \f$r\f$ and \f$g\f$.
 * 
 *  In any case the first thing we need to do is establish the set of real
 *  fragments.  The user does this by setting up another fragmenter and
 *  providing this fragmenter that key via the SYSTEM_FRAGMENTER_KEY option.
 *  For example to do a three-body VMFC computation we set up a fragmenter
 *  to run a three-body fragmentation and provide this module with that key.
 * 
 *  Once we have the real fragments the trick is to now apply the ghost
 *  fragments correctly.  We do this by reading in a map that for a
 *  given \f$r\f$ tells us the values of \f$g\f$ to apply.  For CP this would
 *  be \f$(N-1)\f$ to all \f$r\f$-mers.  For VMFC this is the power set of
 *  \f$(n-r)\f$.  And for MGMBE this is the power set of specified value.  We
 *  can only accept a dictionary of int to int.  If we assume [x]->y means
 *  for real \f$x\f$-mers we take all pairwise unions with the power set of 
 *  y (as ghosts), then VMFC is just [r]->(n-r) and MGMBE is just [r]->g where
 *  g can range from 0 to (N-r) and is specified by the user.  This leaves how
 *  do we treat CP?  By realizing that a MGMBE where for each \f$r\f$ we
 *  pick \f$g=(N-r)\f$ recover a MBE completely in the supersystem basis we 
 *  see how to do this.  This is because all the intermediate quantities
 *  necesarilly cancel (we cheat behind the scenes and don't generate them).
 *  Note however that if for all \f$r\f$ this is not the case, e.g. for monomers
 *  we consider \f$g=(N-1)\f$, but for dimers we pick \f$g=3\f$ then we have
 *  \f[
 *    \sum_{I<J}\sum_{K<L<M}E_{IJ}(KLM)-E_{I}(KLM)-E_J(KLM)
 *  \f]
 *  (fragments in parentheses are ghosts) and:
 *  \f[
 *    \sum_{I<J}\sum_{K<L}E_{IJ}(KL)-E_{I}(KL)-E_J(KL)
 *  \f]
 *  and:
 *  \f[
 *    \sum_{I<J}\sum_{K}E_{IJ}(K)-E_{I}(K)-E_J(K)
 *  \f]
 *  as well as two of these (to avoid double counting):
 *  \f[
 *    \sum_{I<J}E_{IJ}-E_{I}-E_J
 *  \f]
 *  Arguably this is not the greatest choice, but I'm not going to stop you.
 * 
 *  The values of g are read in via the option GHOST_TRUNCATION_ORDERS, which
 *  has the form:
 *  
 *  {{1->2},{2->3}} 
 *  
 *  Which can be read as truncate the monomer g at 2, and the dimer g at 3.  If
 *  there are other n-mers (say trimers) they will not be BSSE corrected.  In
 *  other words, any order not specified is assumed to be non-ghosted.
 *  Because it is annoying to have to know ahead of time how many fragments
 *  there are, we establish that for r->g if(g>=N-r) we assume you want the
 *  supersystem basis applied.
 * 
 *  Implementation details:
 * 
 *  Figuring out the coefficients is the hard part.  For a full CP-like
 *  BSSE correction the weights are whatever they were for the fragments.  For
 *  VMFC they are (-1)^{number of ghost frags mod 2}.  For MGMBE it depends on
 *  the values of g in a non-trivial way.  Although, one could derive a closed
 *  form, such a form is only valid when no n-mers are thrown away.  Hence one
 *  needs to determine them on the fly, much like how we do for a normal MBE.
 *  
 *  To do this realize that the BSSE-free energy \f$\mathcal{E}\f$ can be
 *  written in terms of the supersystem energy \f$E\f$ (possibly determined
 *  with a MBE) and some BSSE correction \f$\delta E_{BSSE}\f$:
 *  \f[
 *  \mathcal{E}=E+\delta E_{BSSE}
 *  \f]
 *  If \f$E\f$ is a MBE then the weights of the n-mers there are whatever
 *  they are normally.  Within the correction we have:
 *  \f[
 *  \delta E_{BSSE}=\sum_{i=1}^n\sum_{i-mers}\left[\Delta\mathcal{E}_{i-mer}-
 *                   \Delta E_{i-mer}\right]
 *  \f]
 *  where the terms on the right are some form of BSSE corrected \f$i\f$-body
 *  interaction and its normal \f$i\f$-body counterpart.  For MGMBE, the quantity
 *  in square brackets is further decomposed into an additional sum running over
 *  the ghost monomers.  Assuming \f$E\f$ is also computed with a MBE that
 *  is truncated at order \f$n\f$, then the normal part of the BSSE correction
 *  cancels with the MBE used for \f$E\f$.  For MGMBE we need to know the 
 *  total number of terms in the decomposition to know the weight of the
 *  normal part.  Furthermore, unlike CP and VMFC the \f$\Delta\mathcal{E}\f$
 *  terms are not disjoint for MGMBE and care has to be taken there as well.
 *  
 *
 */
class Ghoster: public pulsar::modulebase::SystemFragmenter{
    public:
        ///Import the constructor from base class
        using pulsar::modulebase::SystemFragmenter::SystemFragmenter;
        
        ///Actually fragments the molecule
        virtual NMerSetType Fragmentize_(const pulsar::system::System& mol);
};

///Makes a ghoster consistent with a typical CP correction
class CPGhoster: public pulsar::modulebase::SystemFragmenter{
public:
    ///Import the constructor from base class
        using pulsar::modulebase::SystemFragmenter::SystemFragmenter;
        
        ///Actually fragments the molecule
        virtual NMerSetType Fragmentize_(const pulsar::system::System& mol);
};

///Makes a ghoster consistent with a typical VMFC(n) correction
class VMFCGhoster: public pulsar::modulebase::SystemFragmenter{
public:
    ///Import the constructor from base class
        using pulsar::modulebase::SystemFragmenter::SystemFragmenter;
        
        ///Actually fragments the molecule
        virtual NMerSetType Fragmentize_(const pulsar::system::System& mol);
};


#endif /* PULSAR_GHUARD_GHOSTER_HPP */


