import pulsar as psr
import numpy as np
from pulsar_modules import *

class CorrelationEnergy(psr.EnergyMethod):
    """The Helgaker CBS is usually only applied to the correlation energy.
       Given the CORRELATED_KEY and the REFERENCE_KEY takes their difference
       and returns it
    """
    
    def __init__(self,myid):
        super(CorrelationEnergy,self).__init__(myid)
    
    def deriv_(self,order,wfn):
        mod=self.create_child_from_option("CORRELATED_KEY")
        highwfn,highegy=mod.deriv(order,wfn)
        mod=self.create_child_from_option("REFERENCE_KEY")
        lowwfn,lowegy=mod.deriv(order,wfn)
        diff=np.array(highegy)-np.array(lowegy)
        return wfn,diff.tolist()

class HelgakerCBS(psr.EnergyMethod):
  """This module implements the two-point Helgaker extrapolation.

  This module computes arbitrary order derivatives of an energy that has
  been extrapolated via the two-point Helgaker extrapolation technique. Using
  this technique the energy in a basis set of cardinal number \f$X\f$, is:

  \f[
     E_{X}=E_{CBS}+\alpha X^3.
  \f]
  
  For two basis sets of cardinal numbers \f$X\f$ and \f$X+1\f$ we can remove
  \f$\alpha\f$ and solve for the CBS energy obtaining:

  \f[
      E_{CBS}=\frac{E_X(X+1)^3-E_{X+1}X^3}{(X+1)^3-X^3}.
  \f]

  The \f$n\f$-th order derivative of this is simply:

  \f[
    E^{(n)}_{CBS}=\frac{E^{(n)}_X(X+1)^3-E^{(n)}_{X+1}X^3}{(X+1)^3-X^3}
  \f]

  Module Options:
      BASIS_CARDINAL_NUMS (list of int) : The cardinal numbers of the basis 
          sets.

      METHODS (list of str): The keys for the electronic structure method to use                             The first method should be the one associated with
                             the first cardinal number.
  """

  def __init__(self, myid):
    """Registers this module with the ModuleManager.  Internal use only!!!"""
    super(HelgakerCBS, self).__init__(myid)

  def deriv_(self,order,wfn):
      """ The function that computes the derivative

      Args:
         order (int): The order of the derivative to compute.

      Returns:
         The requested derivative as a 1-D list of floats.

      Raises:
          Nothing

      """
      da_methods=self.options().get("METHODS")
      BasisCards=self.options().get("BASIS_CARDINAL_NUMS")
      l13=BasisCards[0]**3;l23=BasisCards[1]**3
      denom=l23-l13
      Cs=[l23/denom,-l13/denom]
      Results=Methods.Methods.run_series_of_methods(
            self.module_manager(),self.id(),da_methods,[wfn],order)
      da_sum=Cs[0]*np.array(Results[0][1])+Cs[1]*np.array(Results[1][1])
      #TODO: Combine wfns
      return wfn,da_sum.tolist()

     
         
