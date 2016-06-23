from .SCF import SCF
from .MP2 import MP2
from .CC import CCSDT
from .CompositeMethods import CP
from .CompositeMethods import HelgakerCBS
from .CompositeMethods import FellerCBS
from .CompositeMethods import FPA
from pulsar.modulemanager import ModuleCreationFuncs


def insert_supermodule():
    cf = ModuleCreationFuncs()
    cf.add_py_creator("SCF", SCF.SCF)
    cf.add_py_creator("MP2",MP2.MP2)
    cf.add_py_creator("CCSD(T)",CCSDT.CCSDT)
    cf.add_py_creator("CP",CP.CP)
    cf.add_py_creator("HelgakerCBS",HelgakerCBS.HelgakerCBS)
    cf.add_py_creator("FellerCBS",FellerCBS.FellerCBS)
    cf.add_py_creator("FPA",FPA.FPA)
    return cf
