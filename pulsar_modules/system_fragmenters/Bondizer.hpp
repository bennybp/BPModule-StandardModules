#pragma once

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief Breaks a system apart into groups of atoms that are seperated by at
 *         most a pre-defined number of bonds
 * 
 *  A detailed description of this class is available [here](@ref bondizer)
 */
class Bondizer : public pulsar::SystemFragmenter
{
public:
    using pulsar::SystemFragmenter::SystemFragmenter;///<Use base constructors

    ///Implements the guts of the SystemFragmenter to be a Bondizer
    virtual pulsar::NMerSetType fragmentize_(const pulsar::System & mol);

};

