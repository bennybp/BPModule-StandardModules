from . import TestPyModule1
from . import TestOptions
from bpmodule.modulemanager import ModuleCreationFuncs


def InsertSupermodule():
    cf = ModuleCreationFuncs()
    cf.AddPyCreator("TestPyModule1", TestPyModule1.TestPyModule1)
    cf.AddPyCreator("TestOptions", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_int", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_float", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_bool", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_str", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_listint", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_listfloat", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_listbool", TestOptions.TestOptions)
    cf.AddPyCreator("TestOptions_liststr", TestOptions.TestOptions)
    return cf
