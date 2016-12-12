import pulsar as psr

def run(mm):
    tester = psr.PyTester("Testing Atomizer")

    water6=psr.make_system("""
    0 1
    O    -2.0449536949999998    -1.6898322539999999     0.0354707364500000
    H    -2.3427132454308994    -2.1474611062791298     0.8216939386571565
    H    -1.1344686596866658    -1.9649570182333860    -0.0720244010028244
    O     2.4859145229999999    -0.9260656876000000     0.0354704040100000
    H     2.2689370278862486    -0.0000001286725659    -0.0720246525395077
    H     3.0311125996609807    -0.9551186438629339     0.8216935421441762
    O    -2.4859145329999999     0.9260657306000000    -0.0354704090300000
    H    -3.0311126049484312     0.9551186852913661    -0.8216935504901171
    H    -2.2689370387823162     0.0000001718485472     0.0720246508442851
    O     2.0449536849999999     1.6898322980000000    -0.0354707319800000
    H     1.1344686495378762     1.9649570616208789     0.0720244057802314
    H     2.3427132361387448     2.1474611537792310    -0.8216939318818128
    O    -0.4409608357000000     2.6158980070000002     0.0354706148100000
    H    -0.6883991903307550     3.1025798186648132     0.8216938113975091
    H    -1.1344683993833504     1.9649572116176803    -0.0720244084283669
    O     0.4409608257000000    -2.6158979640000002    -0.0354706097500000
    H     0.6883991785965635    -3.1025797787366964    -0.8216938049817561
    H     1.1344683896564345    -1.9649571682424622     0.0720244094544075
    """)
    
    mm.load_module("pulsar_modules","Atomizer","PSR_ATOM_FRAG")
    my_mod=mm.get_module("PSR_ATOM_FRAG",0)

    #This will make the right answer
    corr={}
    for i,atomi in enumerate(water6):
        SN=str(i)+" "
        Temp=psr.NMerInfo()
        Temp.nmer=psr.System(water6,False)
        Temp.nmer.insert(atomi)
        Temp.weight=1.0
        Temp.sn=set({i,})
        corr[SN]=Temp
    print("done with corr answer")
    tester.test_return("Resulting fragments are correct",True,
        corr,my_mod.fragmentize,water6)

    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
