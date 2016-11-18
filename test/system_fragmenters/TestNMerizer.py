from TestFxns import *
from itertools import combinations

def make_nmers(weights,corr,waters,n):
    for i,j in corr.items():
        j.weight=weights[len(j.sn)-1]
    for comb in combinations(range(0,6),n):
        key=""
        for i in comb:
            key=key+str(i)+" "
        corr[key]=psr.NMerInfo()
        corr[key].nmer=psr.System(waters[0].nmer,False)
        for i in comb:
            corr[key].nmer.union_assign(waters[i].nmer)
        corr[key].weight=1.0
        corr[key].sn=set(comb)


def run(mm):
    tester = Tester("Testing NMerizer")

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
    
    #This makes the six waters
    waters=[psr.NMerInfo() for x in range(0,6)]
    for i,wateri in enumerate(waters):
        wateri.nmer=psr.System(water6,False)
        wateri.weight=1.0
        wateri.sn={i}
    for i,atomi in enumerate(water6):
        waters[int((i-i%3)/3)].nmer.insert(atomi)
    


    mm.load_module("pulsar_modules","NMerizer","PSR_NMER_FRAG")
    mm.load_module("pulsar_modules","Bondizer","PSR_BOND_FRAG")
    mm.change_option("PSR_NMER_FRAG","SYSTEM_FRAGMENTER_KEY","PSR_BOND_FRAG")
    mm.change_option("PSR_NMER_FRAG","TRUNCATION_ORDER",0)
    
    my_mod=mm.get_module("PSR_NMER_FRAG",0)
    frags=my_mod.fragmentize(water6)
    
    corr={}
    tester.test_value("Resulting fragments are correct n=0",corr,frags)
    my_mod.options().change("TRUNCATION_ORDER",1)
    
    corr={str(i)+" ":waters[i] for i in range(0,6)}
    frags=my_mod.fragmentize(water6)
    tester.test_value("Resulting fragments are correct n=1",corr,frags)
    
    #Weights for n=2 to 6 respectively
    cs=[[-4.0],[6.0,-3.0],[-4.0,3.0,-2.0],[1.0,-1.0,1.0,-1.0],
        [0.0 for i in range(5)]]
    
    for n,c in enumerate(cs):
        make_nmers(c,corr,waters,n+2)
        my_mod.options().change("TRUNCATION_ORDER",n+2)        
        frags=my_mod.fragmentize(water6)
        tester.test_value("Resulting fragments are correct n="+str(n+2),corr,
             frags)
    
    tester.print_results()


with psr.ModuleAdministrator() as mm:
    run(mm)
    
psr.finalize()

