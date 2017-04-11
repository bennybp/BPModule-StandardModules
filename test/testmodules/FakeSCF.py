import pulsar as psr
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
waters=[psr.System(water3,False) for i in range(6)]
for i in range(3):
    for j in range(i*3,i*3+3):waters[i].insert(asu[j])
waters[3]=waters[0]+waters[1]
waters[4]=waters[0]+waters[2]
waters[5]=waters[1]+waters[2]
#Energies are SCF/aug-cc-pvdz
#Respectively monomers 1, 2, 3 dimers 12, 13, 23
Egys=[-76.0414077582648957,-76.0414077580937260,-76.0414077563411723,
      -152.0843137146730442,-152.0850458044913580,-152.0836453173948826]
class FakeSCF(psr.EnergyMethod):
    """Simulates an SCF on water molecules"""
    def __init__(self, myid):
        super(FakeSCF, self).__init__(myid)
    def deriv_(self,order,wfn):
        for i in range(6):
            if wfn.system==waters[i]: return wfn,[Egys[i]]





