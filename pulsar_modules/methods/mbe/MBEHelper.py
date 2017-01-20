import pulsar as psr

def MBE_wrapper(mm,wfn,truncation_order,MBEMod,nmerizer_key,sub_keys=[]):
    """ Relies on Pulsar's cacheing ability to compute 1-body, 2-body, etc.
        derivatives.

        EnergyMethod modules return the best energy they computed.  In the case
        of the MBE this is the n-body energy; however, all the information for
        computing the (n-1)-body, the (n-2)-body, etc. energies is available.
        Within the Pulsar framework the way to get at these energies is to
        simply loop over n requesting the energy (this is effecient because of
        Pulsar's cacheing system).  This function is a wrapper that does just
        that.  Additionally, this function will also, optionally, loop over
        EnergyMethod modules.  This is intended to be used when the target
        computation was something like CCSD(T).   In this case, we know that the
        information required to compute the HF, MP2, and CCSD energies, via the
        MBE, are is available.  Again, within Pulsar, the way to do this is to
        loop over the methods and rely on the cacheing system to do this
        effeciently.

        Args:
        mm: the ModuleManager instance to use
        wfn: the Wavefunction instance to use
        truncation_order: Consider unions of up to how many monomers?
        MBEmod: the MBE EnergyMethod module instance to call
        nmerizer_key: the key for Fragmentizer that makes the n-mers
        sub_keys: a list of EnergyMethod keys to call
    """
    Derivs={meth:{} for meth in sub_keys}
    for i in range(truncation_order,0,-1):
        newkey=nmerizer_key+str(i)
        mm.duplicate_key(nmerizer_key,newkey)
        mm.change_option(newkey,"TRUNCATION_ORDER",i)
        MBEMod.options().change("SYSTEM_FRAGMENTER_KEY",newkey)
        if sub_keys != []:
            for meth in sub_keys:
                MBEMod.options().change("METHOD_KEY",meth)
                NewWfn,Deriv=MBEMod.deriv(0,wfn)
                Derivs[meth][i]=Deriv[0]
        else:
            NewWfn,Deriv=MBEMod.deriv(0,wfn)
            Derivs[i]=Deriv[0]
    return Derivs

