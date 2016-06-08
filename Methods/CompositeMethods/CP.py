import pulsar as psr

class CP(psr.modulebase.EnergyMethod):
   """This module performs a typical CP correction"""

   def __init__(self,myid):
      """Registers this module with ModuleManager.  Internal use only!!!""" 
      super(CP,self).__init__(myid)
   
   def Deriv_(self,order,wfn):
      MethodKey=self.Options().Get("METHOD")
      MBEKey=self.Options().Get("MBE_KEY")
      GhosterKey=self.Options().Get("GHOSTER_KEY")
      NewMBEKey=self.MManager().GenerateUniqueKey()
      self.MManager().DuplicateKey(MBEKey,NewMBEKey)
      self.MManager().ChangeOption(NewMBEKey,"FRAGMENTIZER",GhosterKey)
      
      MyMod=self.CreateChild(MethodKey)
      SSWfn,SSDeriv=MyMod.Deriv(order,wfn)
      MyMod=self.CreateChild(MBEKey)
      MBEWfn,MBEDeriv=MyMod.Deriv(order,wfn)
      MyMod=self.CreateChild(NewMBEKey)
      CPWfn,CPDeriv=MyMod.Deriv(order,wfn)
   
      FinalDeriv=[SSDeriv[i]-CPDeriv[i]+MBEDeriv[i] for i in range(0,len(SSDeriv))]
      return wfn,FinalDeriv