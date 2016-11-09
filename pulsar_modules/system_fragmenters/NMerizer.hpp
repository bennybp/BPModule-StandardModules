#pragma once

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief Takes unions of up to \f$N\f$ fragments
 * 
 *  A detailed description of this class is available [here](@ref nmerizer)
 */
class NMerizer : public pulsar::SystemFragmenter
{
public:
    using pulsar::SystemFragmenter::SystemFragmenter;///<Use base constructors

    ///Implements the guts of the SystemFragmenter to be an NMerizer
    virtual pulsar::NMerSetType fragmentize_(const pulsar::System & mol);

};
