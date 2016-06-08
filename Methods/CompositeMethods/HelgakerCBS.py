import pulsar as psr

class HelgakerCBS(psr.modulebase.EnergyMethod):
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

  The @f$n@f$-th order derivative of this is simply:

  \f[
    E^{(n)}_{CBS}=\frac{E^{(n)}_X(X+1)^3-E^{(n)}_{X+1}X^3}{(X+1)^3-X^3}
  \f]

  Module Options:
      BASIS_CARDINAL_NUMS (list of int) : The cardinal numbers of the basis 
          sets.

      METHOD (str): The key for the electronic structure method to use.
  """

  def __init__(self, myid):
    """Registers this module with the ModuleManager.  Internal use only!!!"""
    super(HelgakerCBS, self).__init__(myid)

  def Deriv_(self,order,wfn):
      """ The function that computes the derivative

      Args:
         order (int): The order of the derivative to compute.

      Returns:
         The requested derivative as a 1-D list of floats.

      Raises:
          Nothing

      """
      BasisCards=self.Options().Get("BASIS_CARDINAL_NUMS")
      l13=BasisCards[0]**3
      l23=BasisCards[1]**3
      denom=l23-l13
      Cs=[l23/denom,-l13/denom]
      MIM=self.CreateChild(self.Options().Get("MIM_KEY"))
      MIM.Options().Change("WEIGHTS",Cs)
      return MIM.Deriv(order,wfn)

     
         
