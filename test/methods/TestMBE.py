import pulsar as psr

class TestMethod(psr.EnergyMethod):
  def __init__(self, myid):
    super(TestMethod, self).__init__(myid)

  def deriv_(self,order,wfn):
      deriv=[1.0]
      if order==0:
          for ai in wfn.system:
              for qi in ai:
                  deriv[0]*=qi
      return wfn,deriv


def run(mm,tester):
    water=psr.make_system("""
    0 1
    O    -2.0449536949999998    -1.6898322539999999     0.0354707364500000
    H    -2.3427132454308994    -2.1474611062791298     0.8216939386571565
    H    -1.1344686596866658    -1.9649570182333860    -0.0720244010028244
    """)
    mm.load_module("pulsar_modules","NMerizer","PSR_NMER_FRAG")
    mm.load_module("pulsar_modules","Atomizer","PSR_BOND_FRAG")
    mm.load_module("pulsar_modules","MBE","PSR_MBE")
    minfo=psr.ModuleInfo()
    minfo.name="FakeEnergyMethod"
    minfo.type="python_module"
    minfo.path="./TestMBE.py"
    minfo.base="EnergyMethod"
    minfo.version="1.0"
    minfo.description="A 3*Natoms H.O."
    mm.load_module_from_minfo(minfo,"PSR_EGY_METH")
    mm.change_option("PSR_NMER_FRAG","SYSTEM_FRAGMENTER_KEY","PSR_BOND_FRAG")
    mm.change_option("PSR_NMER_FRAG","TRUNCATION_ORDER",1)
    mm.change_option("PSR_MBE","SYSTEM_FRAGMENTER_KEY","PSR_NMER_FRAG")
    mm.change_option("PSR_MBE","METHOD_KEY","PSR_EGY_METH")
    mm.change_option("PSR_EGY_METH","MAX_DERIV",1)
    
    my_mod=mm.get_module("PSR_MBE",0)
    
    #Correct answers
    carts=[],egy=[1.0 for ai in water]
    for i,ai in enumerate(water):
        for q in ai:
            carts.append(q)
            egy[i]*=q
    E1b=egy[0]+egy[1]+egy[2]
    dE1b=[]
    for i in range(3):
        dE1b.append(egy[0]/carts[i]+egy[1]+egy[2])
    for i in range(3,6):
        dE1b.append(egy[0]+egy[1]/carts[i]+egy[2])
    for i in range(6,9):
        dE1b.append(egy[0]+egy[1]+egy[2]/carts[i])
            
    E2b=egy[0]*egy[1]+egy[0]*egy[2]+egy[1]*egy[2]
    
    deriv=my_mod.deriv(0,wfn)
    print(deriv[1])
        
        
    
    

with psr.ModuleAdministrator() as mm:
    run(mm,tester)
psr.finalize()
tester.print_results()
