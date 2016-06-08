from .SCF import SCF
from .MP2 import MP2
from .CC import CCSDT
from .CompositeMethods import CP
from .CompositeMethods import HelgakerCBS
from .CompositeMethods import FellerCBS
from .CompositeMethods import FPA
from pulsar.modulemanager import ModuleCreationFuncs


def InsertSupermodule():
    cf = ModuleCreationFuncs()
    cf.AddPyCreator("SCF", SCF.SCF)
    cf.AddPyCreator("MP2",MP2.MP2)
    cf.AddPyCreator("CCSD(T)",CCSDT.CCSDT)
    cf.AddPyCreator("CP",CP.CP)
    cf.AddPyCreator("HelgakerCBS",HelgakerCBS.HelgakerCBS)
    cf.AddPyCreator("FellerCBS",FellerCBS.FellerCBS)
    cf.AddPyCreator("FPA",FPA.FPA)
    return cf
