import os
import sys
import pulsar as psr
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
print(sys.path)
from pulsar_modules.analysis.MBEHelper import *

def run(mm):
    tester=psr.PyTester("Testing the MBEHelper Analysis Module")
    ang2au=1/0.52917721067
    carts=[[-2.0449536949999998,-1.6898322539999999,0.0354707364500000],
           [-2.3427132454308994,-2.1474611062791298,0.8216939386571565],
           [-1.1344686596866658,-1.9649570182333860,-0.0720244010028244],
           [2.4859145229999999,-0.9260656876000000,0.0354704040100000],
           [2.2689370278862486,-0.0000001286725659,-0.0720246525395077],
           [3.0311125996609807,-0.9551186438629339,0.8216935421441762],
           [-2.4859145329999999,0.9260657306000000,-0.0354704090300000],
           [-3.0311126049484312,0.9551186852913661,-0.8216935504901171],
           [-2.2689370387823162,0.0000001718485472,0.0720246508442851]
    ]
    for i in carts:
        for j in range(3):i[j]*=ang2au

    atoms=[
      psr.create_atom(carts[0],8),
      psr.create_atom(carts[1],1),
      psr.create_atom(carts[2],1),
      psr.create_atom(carts[3],8),
      psr.create_atom(carts[4],1),
      psr.create_atom(carts[5],1),
      psr.create_atom(carts[6],8),
      psr.create_atom(carts[7],1),
      psr.create_atom(carts[8],1)
    ]

    asu=psr.AtomSetUniverse()
    for ai in atoms:asu.insert(ai)
    water3=psr.System(asu,True)
    mm.load_module("pulsar_modules","NMerizer","PSR_NMER_FRAG")
    mm.load_module("pulsar_modules","Bondizer","PSR_BOND_FRAG")
    mm.load_module("pulsar_modules","MBE","PSR_MBE")
    mm.load_module("testmodules","Fake SCF","Fake SCF")
    mm.change_option("PSR_NMER_FRAG","SYSTEM_FRAGMENTER_KEY","PSR_BOND_FRAG")
    mm.change_option("PSR_NMER_FRAG","TRUNCATION_ORDER",2)
    mm.change_option("PSR_MBE","SYSTEM_FRAGMENTER_KEY","PSR_BOND_FRAG")
    mm.change_option("PSR_MBE","METHOD_KEY","Fake SCF")

    my_mod=mm.get_module("PSR_MBE",0)
    wfn=psr.Wavefunction()
    wfn.system=water3
    egy=mbe_wrapper(mm,wfn,2,my_mod,"PSR_NMER_FRAG",['Fake SCF'])
    corr_egys={1:-228.1242232726998,2:-228.1287815638595}
    for meth,test_egy in egy.items():
        for order,egy_i in test_egy.items():
            msg=meth+" "+str(order)+" -body energy"
            tester.test_double(msg,egy_i,corr_egys[order])

    corr_ints={'1 ': [-76.04140775809373], '0 1 ': [-0.0014981983144224387],
               '0 2 ': [-0.002230289885289949], '0 ': [-76.0414077582649],
               '1 2 ': [-0.0008298029599842494], '2 ': [-76.04140775634117]}
    ints=mbe_interactions(mm,wfn,2,my_mod,"PSR_NMER_FRAG",['Fake SCF'])
    for meth,test_int in ints.items():
        for inti,egyi in test_int.items():
            msg=meth+" "+str(inti)+" "+str(len(inti.split()))+"-body interaction"
            tester.test_double(msg,egyi,corr_ints[inti])

    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
