import pulsar as psr
import itertools
import numpy as np

class MyCzyCM(psr.EnergyMethod):
    """This is an example of how to nest a bunch of methods
    
       This method will perform a focal point analysis where the HF energy is
       computed in the aQZ basis set, the MP2 correlation energy is the result
       of an aTZ,aQZ CBS extrapolation and the CCSD(T) correlation energy is
       the result of the CBS extrapolation of an aDZ, aTZ extrapolation
       
       That is:
       HF/aQZ+MP2/a(TQ)Z+\f$\delta\f$CCSD(T)/a(DT)Z
    """
    
    def __init__(self, myid):
        """Registers this module with the ModuleManager.  
           Internal use only!!!
        """
        super(MyCzyCM, self).__init__(myid)
    def deriv_(self,order,wfn):
        mm=self.module_manager()
        #First set up the CCSD(T)
        meths=["CCSD(T)","MP2","SCF","MP2_DRY"]
        bs=["aug-cc-pvdz","aug-cc-pvtz","aug-cc-pvqz"]
        for lot in itertools.product(meths,bs):
            key=lot[0]+lot[1]
            mm.duplicate_key("PSI4_"+lot[0],key)
            mm.change_option(key,"BASIS_SET",lot[1])
            mm.duplicate_key("PSR_CORRELATION_ENERGY",key+"CE")
            mm.change_option(key+"CE","CORRELATED_KEY",key)
            mm.change_option(key+"CE","REFERENCE_KEY","SCF"+lot[1])
        mm.duplicate_key("PSR_HELGAKER_CBS","CCSD(T)/CBS")
        mm.duplicate_key("PSR_HELGAKER_CBS","MP2/a(DT)Z")
        mm.duplicate_key("PSR_HELGAKER_CBS","MP2/a(TQ)Z")
        
        psr.datastore.set_options(mm,
            {"CCSD(T)/CBS":{
                "METHODS":["CCSD(T)aug-cc-pvdzCE","CCSD(T)aug-cc-pvtzCE"],
                "BASIS_CARDINAL_NUMS":[2,3]},
             "MP2/a(DT)Z":{
                "METHODS":["MP2aug-cc-pvdzCE","MP2aug-cc-pvtzCE"],
                "BASIS_CARDINAL_NUMS":[2,3]},
             "MP2/a(TQ)Z":{
                "METHODS":["MP2aug-cc-pvtzCE","MP2aug-cc-pvqzCE"],
                "BASIS_CARDINAL_NUMS":[3,4]},
             "PSR_FPA":{"LARGE_MP2_KEY":"MP2/a(TQ)Z",
                        "SMALL_MP2_KEY":"MP2/a(DT)Z",
                        "CCSD(T)_KEY":"CCSD(T)/CBS"
             }
            })
        mod=self.create_child("PSR_FPA")
        ce_wfn,ce_egy=mod.deriv(order,wfn)
        mod=self.create_child("SCFaug-cc-pvqz")
        hf_wfn,hf_egy=mod.deriv(order,wfn)

        da_sum=np.array(hf_egy)+np.array(ce_egy)
        #TODO: Combine Wfns
        return wfn,da_sum.tolist()
    
