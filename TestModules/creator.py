from . import TestPyModule1
from . import TestOptions
from pulsar.modulemanager import ModuleCreationFuncs


def insert_supermodule():
    cf = ModuleCreationFuncs()
    cf.add_py_creator("TestPyModule1", TestPyModule1.TestPyModule1)
    cf.add_py_creator("TestOptions", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_int", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_float", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_bool", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_str", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_listint", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_listfloat", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_listbool", TestOptions.TestOptions)
    cf.add_py_creator("TestOptions_liststr", TestOptions.TestOptions)
    return cf
