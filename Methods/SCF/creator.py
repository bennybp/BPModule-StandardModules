from . import SCF
from bpmodule.modulemanager import ModuleCreationFuncs


def InsertSupermodule():
    cf = ModuleCreationFuncs()
    cf.AddPyCreator("SCF", SCF.SCF)
    return cf