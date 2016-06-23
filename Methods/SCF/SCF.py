import pulsar as psr
from subprocess import call
from CallPsi4 import CallPsi4

class SCF(psr.modulebase.EnergyMethod):
  def __init__(self, myid):
    super(SCF, self).__init__(myid)

  def deriv_(self,order,wfn):
      return CallPsi4(self,"SCF",order,wfn)