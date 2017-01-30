import pulsar as psr
import copy
from itertools import combinations
def mbe_wrapper(mm,wfn,truncation_order,MBEMod,nmerizer_key,sub_keys=[]):
    """ Relies on Pulsar's cacheing ability to compute 1-body, 2-body, etc.
        derivatives.

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
        if not mm.has_key(newkey):
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

def mbe_interactions(mm,wfn,truncation_order,MBEMod,nmerizer_key,sub_keys=[]):
    """Relies on Pulsar's cacheing ability to compute individual 1-body, 2-body,
       etc. interactions.

       Args:
       mm: the ModuleManager instance to use
       wfn: the Wavefunction instance to use
       truncation_order: Consider unions of up to how many monomers?
       MBEmod: the MBE EnergyMethod module instance to call
       nmerizer_key: the key for Fragmentizer that makes the n-mers
       sub_keys: a list of EnergyMethod keys to call
    """
    results={meth:{} for meth in sub_keys}
    newkey=nmerizer_key+str(truncation_order)
    if not mm.has_key(newkey):
        mm.duplicate_key(nmerizer_key,newkey)
        mm.change_option(newkey,"TRUNCATION_ORDER",i)

    my_mod=mm.get_module(newkey,0)
    frags= my_mod.fragmentize(wfn.system)
    newwfn=copy.deepcopy(wfn)
    for sn,mol in frags.items():
        newwfn.system=psr.System(mol.nmer)
        for meth in sub_keys:
            egy_meth=mm.get_module(meth,0)
            results[meth][sn]=egy_meth.deriv(0,newwfn)[1]
    ints=copy.deepcopy(results)
    for meth,egys in results.items():
        for sn,egy in egys.items():
            ei=ints[meth][sn]
            for n in range(1,len(sn.split())):
                for c in combinations(sn.split(),n):
                    key=""
                    for k in c:key=key+k+" "
                    ei[0]+=(-1)**((len(sn)-n)%2)*egys[key][0]
    return ints
