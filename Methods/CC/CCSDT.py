import pulsar as psr
from subprocess import call
from CallPsi4 import CallPsi4

class CCSDT(psr.modulebase.EnergyMethod):
  def __init__(self, myid):
    super(CCSDT, self).__init__(myid)

  def Deriv_(self,order,wfn):
      return CallPsi4(self,"CCSD(T)",order,wfn)
