import bpmodule as bp

class FellerCBS(bp.modulebase.EnergyMethod):
  """This module implements the three-point Feller extrapolation

  In the Feller extrapolation scheme, the energy for a basis set of cardinal 
  number \f$X\f$ is given by:
  \f[
    E_X=E_{CBS}+\alpha\exp(-\beta X)
  \f]
  For three basis sets of cardinal numbers: \f$X\f$,\f$X+1\f$, and \f$X+2\f$,
  the equations can be rearranged to give:
  \f[
     E_{CBS}=E_X-\frac{\left(E_X-E_{X+1}\right)^2}{E_X-2E_{X+1}-E_{X+2}}
  \f]

  """
  def __init__(self, myid):
    super(HelgakerCBS, self).__init__(myid)

  def Deriv_(self,order):
      BasisCards=Options().Get("BASIS_CARDINAL_NUMS")
      l13=BasisCards[0]**3
      l23=BasisCards[1]**3
      denom=l23-l13
      Cs=[l23/denom,-l13/denom]
      MIM=self.GetSubModule("MIM")
      MIM.Options().Change("WEIGHTS",Cs)
      MIM.Options().Change("METHODS",[Options.Get("METHOD")])
      return MIM.Deriv(order)
