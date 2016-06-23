import pulsar as psr

class CP(psr.modulebase.EnergyMethod):
   """This module performs a typical CP correction"""

   def __init__(self,myid):
      """Registers this module with ModuleManager.  Internal use only!!!""" 
      super(CP,self).__init__(myid)
   
   def deriv_(self,order,wfn):
      MethodKey=self.options().get("METHOD")
      MBEKey=self.options().get("MBE_KEY")
      GhosterKey=self.options().get("GHOSTER_KEY")
      NewMBEKey=self.module_manager().generate_unique_key()
      self.module_manager().duplicate_key(MBEKey,NewMBEKey)
      self.module_manager().change_option(NewMBEKey,"FRAGMENTIZER",GhosterKey)
      
      MyMod=self.create_child(MethodKey)
      SSWfn,SSDeriv=MyMod.deriv(order,wfn)
      MyMod=self.create_child(MBEKey)
      MBEWfn,MBEDeriv=MyMod.deriv(order,wfn)
      MyMod=self.create_child(NewMBEKey)
      CPWfn,CPDeriv=MyMod.deriv(order,wfn)
   
      FinalDeriv=[SSDeriv[i]-CPDeriv[i]+MBEDeriv[i] for i in range(0,len(SSDeriv))]
      return wfn,FinalDeriv