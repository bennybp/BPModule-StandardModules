from .FakeSCF import *
from pulsar import ModuleCreationFuncs

def insert_supermodule():
    cf = ModuleCreationFuncs()
    cf.add_py_creator("Fake SCF", FakeSCF)
    return cf
