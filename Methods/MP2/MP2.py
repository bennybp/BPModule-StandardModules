import bpmodule as bp
from subprocess import call
from CallPsi4 import CallPsi4

class MP2(bp.modulebase.EnergyMethod):
  def __init__(self, myid):
    super(MP2, self).__init__(myid)

  def Deriv_(self,order):
      return CallPsi4(self,"MP2",order)
           

