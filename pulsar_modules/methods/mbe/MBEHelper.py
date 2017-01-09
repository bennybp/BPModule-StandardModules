import pulsar as psr

def MBE_wrapper(mm,wfn,truncation_order,MBEMod,nmerizer_key,sub_keys=[]):
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

