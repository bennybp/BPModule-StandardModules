/** \file CrystalFragger.hpp
 *  \brief Write Me!!!!!!
 *  \author Ryan M. Richard
 *  \version 1.0
 *  \date May 26, 2016
 */

#pragma once

#include <pulsar/modulebase/SystemFragmenter.hpp>

/** \brief This class is in charge of fragmenting a periodic system 
 * 
 * 
 *  Owing to the periodic nature of crystals we have to do things a bit
 *  different than a normal system fragmenter. The input to this fragmenter
 *  is a unit cell.  From the option LATTICE_CUTOFF you tell us how far out
 *  you want your n-tuples to be considered (from the unit cell), thus a selection
 *  of [1,1,1] would mean n-tuples are taken from all unit cells at most 1 cell
 *  away (i.e. we are considering a 3 x 3 x 3 supercell).  We restrict each
 *  n-tuple to having at least one monomer in the unit cell [note that the
 *  the intersections of these n-tuples, the (n-1)-tuples, the (n-2)-tuples,
 *  etc. necessarilly do not satisfy this criteria (consider a trimer that
 *  has one monomer in the UC, two of its constitieunt dimers also have one
 *  monomer in the UC, but the other does not.  Neglecting this latter dimer
 *  will lead to an incorrect three-body energy and can't be done).
 * 
 *  Once we have generated all n-mers and their intersections, the challenge 
 *  becomes determining the unique ones and setting the coefficients correctly.
 *  To do this we steal Beran's strategy and superimpose the n-mers, if they
 *  are the same we only keep one and up its coefficient.
 * 
 */
class CrystalFragger : public pulsar::SystemFragmenter{
public:
    ///Import the constructor from the base class
    using pulsar::SystemFragmenter::SystemFragmenter;
    
    virtual pulsar::NMerSetType fragmentize_(const pulsar::System & mol);
};

