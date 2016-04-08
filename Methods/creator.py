from .SCF import SCF
from .CBS import Helgaker
from bpmodule.modulemanager import ModuleCreationFuncs


def InsertSupermodule():
    cf = ModuleCreationFuncs()
    cf.AddPyCreator("SCF", SCF.SCF)
    cf.AddPyCreator("HelgakerCBS",Helgaker.HelgakerCBS)
    return cf
