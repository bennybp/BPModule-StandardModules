from .SCF import SCF
from .MP2 import MP2
from .CC import CCSDT
from .CBS import Helgaker
from .CBS import Feller
from .MIM import FPA
from bpmodule.modulemanager import ModuleCreationFuncs


def InsertSupermodule():
    cf = ModuleCreationFuncs()
    cf.AddPyCreator("SCF", SCF.SCF)
    cf.AddPyCreator("MP2",MP2.MP2)
    cf.AddPyCreator("CCSD(T)",CCSDT.CCSDT)
    cf.AddPyCreator("HelgakerCBS",Helgaker.HelgakerCBS)
    cf.AddPyCreator("FellerCBS",Feller.FellerCBS)
    cf.AddPyCreator("FPA",FPA.FPA)
    return cf
