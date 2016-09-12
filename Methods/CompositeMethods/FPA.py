import pulsar as psr
import Methods.Methods
import numpy as np

class FPA(psr.modulebase.EnergyMethod):
   """This module performs a focal point analysis (FPA)
   
      For the purposes of this module an FPA is defined as:
      \f[
         E^{large}=E_{MP2}^{large}+E_{CCSD(T)}^{small}-E_{MP2}^{small},
      \f]
      where "large" and "small" refer to the relative sizes of the basis sets
      used.  Technically this is one particular flavor of FPA, but it is one
      of the more common ones.

      Module Options:
        LARGE_MP2_KEY (str) : What module key should be used to find the MP2
                              module for the large basis set?

        SMALL_MP2_KEY (str) : What module key should be used to find the MP2
                              module for the small basis set?  Defaults to
                              "BP_MP2"
       
        CCSD(T)_KEY (str) : What module key should be used to compute the
                            CCSD(T) energy?  Defaults to "BP_CCSD(T)".  For
                            most cases the MP2_KEY of whatever module 
                            "BP_CSSD(T)" maps to should be the same as 
                            SMALL_MP2_KEY, since the MP2 is free from CCSD(T)
   """

   def __init__(self,myid):
      """Registers this module with ModuleManager.  Internal use only!!!""" 
      super(FPA,self).__init__(myid)
   def deriv_(self,order,wfn):
      """ The function that computes the derivative

      Args:
         order (int): The order of the derivative to compute.

      Returns:
         The requested derivative as a 1-D list of floats.

      Raises:
          Nothing

      """
      Results=Methods.Methods.run_series_of_methods(
                    self.module_manager(),self.id(),
                    [self.options().get("CCSD(T)_KEY"),
                     self.options().get("LARGE_MP2_KEY"),
                     self.options().get("SMALL_MP2_KEY")],
                     [wfn],
                     order)
      da_sum=np.array(Results[0][1])+np.array(Results[1][1])-np.array(Results[2][1])
      #TODO: Combine wfns
      return wfn,da_sum.tolist()
