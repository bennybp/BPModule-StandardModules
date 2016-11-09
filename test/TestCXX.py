from TestFxns import *
import sys
import os

Test_Name=sys.argv[1]
mod_path=os.path.dirname(os.path.realpath(__file__))
mod_path=os.path.join(mod_path,Test_Name+".so")

def run():
    minfo = psr.ModuleInfo()
    minfo.name=Test_Name
    minfo.base="Test_Base"
    minfo.path=mod_path
    minfo.type="c_module"
    minfo.version="0.1a"
    minfo.description="Tests Stuff"
    mm.load_module_from_minfo(minfo,"TestPoint")
    my_mod=mm.get_module("TestPoint",0)
    my_mod.run_test()


with psr.ModuleAdministrator() as mm:
    run()
psr.finalize()

