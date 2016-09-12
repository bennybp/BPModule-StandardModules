from .CompositeMethods import CP
from .CompositeMethods import HelgakerCBS
from .CompositeMethods import FellerCBS
from .CompositeMethods import FPA
from .CompositeMethods import MyCrazyCompositeMethod
from pulsar.modulemanager import ModuleCreationFuncs


def insert_supermodule():
    cf = ModuleCreationFuncs()
    cf.add_py_creator("CP",CP.CP)
    cf.add_py_creator("CorrelationEnergy",HelgakerCBS.CorrelationEnergy)
    cf.add_py_creator("HelgakerCBS",HelgakerCBS.HelgakerCBS)
    cf.add_py_creator("FellerCBS",FellerCBS.FellerCBS)
    cf.add_py_creator("FPA",FPA.FPA)
    cf.add_py_creator("MyCrzyCompMeth",MyCrazyCompositeMethod.MyCzyCM)
    return cf
