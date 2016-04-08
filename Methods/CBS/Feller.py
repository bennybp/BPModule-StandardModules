import bpmodule as bp

class FellerCBS(bp.modulebase.EnergyMethod):
  """This module implements the three-point Helgaker
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
