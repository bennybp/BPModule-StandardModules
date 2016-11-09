from TestFxns import *

def run(mm):
    tester = Tester("Testing Bondizer")

    water2=psr.make_system("""
    0 1
    O    -0.0638166135700000     1.3782116490000000     0.0000013991477170
    H     0.8583313934960812     1.6348726952004131    -0.0000192736448436
    H    -0.0465523195463249     0.4211673531851549     0.0000010357024237
    O     0.0596797744600000    -1.2953846140000000    -0.0000013030828870
    H    -0.3727778924225715    -1.6906360791564299     0.7569587260112379
    H    -0.3728118210906938    -1.6906361161127510    -0.7569419281358893
    """)
    
    mm.load_module("pulsar_modules","Bondizer","PSR_BOND_FRAG")
    my_mod=mm.get_module("PSR_BOND_FRAG",0)
    frags=my_mod.fragmentize(water2)

    #This will make the right answer
    Frag1,Frag2=psr.NMerInfo(),psr.NMerInfo()
    Frag1.nmer,Frag2.nmer=psr.System(water2,False),psr.System(water2,False)
    for i,ai in enumerate(water2):
        if i<3: Frag1.nmer.insert(ai)
        else: Frag2.nmer.insert(ai)
    Frag1.weight,Frag2.weight=1.0,1.0
    Frag1.sn,Frag2.sn={0},{1}
    corr={str(0)+" ":Frag1,str(1)+" ":Frag2}  
    
    tester.test_value("Resulting water dimer fragments are correct",corr,frags)
    
    ethane=psr.make_system("""
    C 0.00 0.00 0.00
    C 0.00 0.00 1.52
    H 1.02 0.00 -0.39
    H -0.51 -0.88 -0.39
    H -0.51 0.88 -0.39
    H -1.02 0.00 1.92
    H 0.51 -0.88 1.92
    H 0.51 0.88 1.92
    """)
    my_mod.options().change("MAX_NBONDS",2)
    frags=my_mod.fragmentize(ethane)
    Frag3=psr.NMerInfo()
    Frag1.nmer,Frag2.nmer,Frag3.nmer=[psr.System(ethane,False) for i in range(0,3)]
    for i,ai in enumerate(ethane):
        Frag1.nmer.insert(ai)
        if i in [0,1,2,3,4]:
            Frag2.nmer.insert(ai)
        if i in [0,1,5,6,7]:
            Frag3.nmer.insert(ai)
    Frag3.weight=1.0
    Frag3.sn={2}
    corr["0 "]=Frag1
    corr["1 "]=Frag2
    corr["2 "]=Frag3
       
    tester.test_value("Ethane fragmented correctly",corr,frags)

    tester.print_results()


with psr.ModuleAdministrator() as mm:
    run(mm)
    
psr.finalize()