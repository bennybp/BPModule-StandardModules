import pulsar as psr

class EEQMMM(psr.EnergyMethod):
    """This module performs an electrostatic embedding QM/MM computation.

       This module assumes that the QM and MM regions are disjoint and that
       they are not bonded.
    """

    def __init__(self,myid):
       """Registers this module with ModuleManager.  Internal use only!!!"""
       super(EEQMMM,self).__init__(myid)

    def deriv_(self,order,wfn):
        qm_mod=self.create_child_from_option("QM_METHOD_KEY")
        mm_mod=self.create_child_from_option("MM_METHOD_KEY")
        mm_chg=self.create_child_from_option("MM_CHARGES_KEY")
        partition=self.create_child_from_option("SYSTEM_FRAGMENTER_KEY")
        qm_key=str(self.options().get("QM_REGION_KEY"))+' '
        mm_key=str(self.options().get("MM_REGION_KEY"))+' '
        systems=partition.fragmentize(wfn.system)
        qm_u = systems[qm_key].nmer.as_universe()
        mm_sys = psr.System(systems[mm_key].nmer)
        wfn.system=mm_sys
        mm_wfn,mm_egy = mm_mod.deriv(order,wfn)
        qs=mm_chg.calculate(0,wfn)
        for q,a in zip(qs,mm_sys):
            qm_u.insert(psr.make_point_charge(a,q))
        wfn.system=psr.System(qm_u,True)
        qm_wfn,qm_egy = qm_mod.deriv(order,wfn)
        wfn.system=psr.System(mm_wfn.system.as_universe()+qm_wfn.system.as_universe(),True)
        return wfn,[mm_egy[i]+qm_egy[i] for i in range(len(qm_egy))]


"""
Psi4 can not handle non-integer nuclear charge, will need to test with another
program
"""
#class JanusQMMM(psr.EnergyMethod):
#   """This module performs a Janus electrostatic embedding QM/MM computation"""
#
#   def __init__(self,myid):
#      """Registers this module with ModuleManager.  Internal use only!!!"""
#      super(JanusQMMM,self).__init__(myid)
#
#   def deriv_(self,order,wfn):
#       qm_mod=self.create_child_from_option("QM_METHOD_KEY")
#       mm_mod=self.create_child_from_option("MM_METHOD_KEY")
#       mm_chg=self.create_child_from_option("MM_CHARGES_KEY")
#       partition=self.create_child_from_option("SYSTEM_FRAGMENTER_KEY")
#       qm_key=str(self.options().get("QM_REGION_KEY"))+' '
#       mm_key=str(self.options().get("MM_REGION_KEY"))+' '
#       yy_key=self.options().get("YING_YANG_KEY")
#       systems=partition.fragmentize(wfn.system)
#       qm_u = systems[qm_key].nmer.as_universe()
#       wfn.system=psr.System(systems[mm_key].nmer,True)
#       #mm_wfn,mm_egy = mm_mod.deriv(order,wfn)
#       #qs=mm_chg.calculate(0,wfn)
#       #print(mm_egy)
#       #print(qs)
#       wfn.system=psr.System(qm_u,True)
#       qm_wfn,qm_egy = qm_mod.deriv(order,wfn)
#       #wfn.system=psr.System(mm_wfn.system.as_universe()+qm_wfn.system.as_universe(),True)
#       return wfn,[0.0]

