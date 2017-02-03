import os
import sys
import math
import pulsar as psr
from pulsar_modules.analysis.PESScan import make_pes_range,pes_scan

corr_dist_scan=[
psr.make_system("""
O             -0.06381661        1.37821165        0.00000140
H              0.85833139        1.63487270       -0.00001927
H             -0.04655232        0.42116735        0.00000104
O             -0.20328150       -2.26019100        0.00000461
H             -0.63573916       -2.65544247        0.75696464
H             -0.63577309       -2.65544251       -0.75693601
"""),
psr.make_system("""
O             -0.06381661        1.37821165        0.00000140
H              0.85833139        1.63487270       -0.00001927
H             -0.04655232        0.42116735        0.00000104
O             -0.46624277       -3.22499739        0.00001053
H             -0.89870043       -3.62024886        0.75697056
H             -0.89873436       -3.62024890       -0.75693009
""")
]

corr_angle_scan=[
psr.make_system("""
H              0.85833139        1.63487270       -0.00001927
O             -0.06381661        1.37821165        0.00000140
H             -0.04655232        0.42116735        0.00000104
"""),
psr.make_system("""
H              0.83245285        1.71426652       -0.00001869
O             -0.06381661        1.37821165        0.00000140
H             -0.04655232        0.42116735        0.00000104
""")
]

corr_torsion_scan=[
psr.make_system("""
H              1.18510000       -0.00390000        0.98750000
C              0.75160000       -0.02250000       -0.02090000
H              1.16690000        0.83300000       -0.56930000
H              1.11550000       -0.93290000       -0.51450000
C             -0.75160000        0.02250000        0.02090000
H             -1.16690000       -0.83340000        0.56870000
H             -1.11570000        0.93260000        0.51510000
H             -1.18500000        0.00440000       -0.98750000
"""),
psr.make_system("""
H              1.17950737       -0.18141875        0.97748814
C              0.75160000       -0.02250000       -0.02090000
H              1.17370825        0.91296476       -0.41055011
H              1.11428371       -0.83536722       -0.66323927
C             -0.75160000        0.02250000        0.02090000
H             -1.16690000       -0.83340000        0.56870000
H             -1.11570000        0.93260000        0.51510000
H             -1.18500000        0.00440000       -0.98750000
""")
]

corr_imp_scan=[
psr.make_system("""
H              0.26580000        0.64960000        0.68220000
N             -0.03530000       -0.04400000        0.02850000
H              0.77740000       -0.45320000       -0.38500000
H             -0.55220000        0.41480000       -0.69350000
"""),
psr.make_system("""
H              0.25404618        0.75492327        0.55487084
N             -0.03530000       -0.04400000        0.02850000
H              0.77740000       -0.45320000       -0.38500000
H             -0.55220000        0.41480000       -0.69350000
""")
]

def run_scan(systems,points,scan_range,corr,desc,tester):
    for i,sys in enumerate(pes_scan(systems,points,scan_range)):
        pes_msg=desc+" PES point "+str(i)
        for j,(atom1,atom2) in enumerate(zip(sys,corr[i])):
            atm_msg=" atom "+str(j)
            for k,(q1,q2) in enumerate(zip(atom1,atom2)):
                tester.test_double(pes_msg+atm_msg+" coord: "+str(k),q1,q2)


def run(mm):
    tester=psr.PyTester("Testing the PESScan Analysis Module")
    water2=corr_dist_scan[0].as_universe()

    ###Distance scan HO-H----OH2 distance
    systems=[psr.System(water2,False) for i in range(2)]
    for i in range(3):
        systems[0].insert(water2[i])
        systems[1].insert(water2[i+3])

    points=[water2[1],water2[3]]
    scan_range=make_pes_range(0,1,1,1/0.52917721067)
    run_scan(systems,points,scan_range,corr_dist_scan,"Distance stretch",tester)


    ###Angle scan of water H-O-H angle
    new_systems=[psr.System(systems[0],False) for i in range(2)]
    points=[water2[1],water2[0],water2[2]]
    new_systems[0].insert(points[0])
    new_systems[1].insert(points[1])
    new_systems[1].insert(points[2])
    scan_range=[i*math.pi/180.0 for i in range(0,6,5)]
    run_scan(new_systems,points,scan_range,corr_angle_scan,"Angle bend",tester)


    ###Torsion scan of ethane H-CC-H angle
    ethane=corr_torsion_scan[0].as_universe()
    systems=[psr.System(ethane,False) for i in range(2)]
    for i in range(4):
        systems[0].insert(ethane[i])
        systems[1].insert(ethane[i+4])
    points=[ethane[0],ethane[1],ethane[4],ethane[5]]
    scan_range=[i*math.pi/180.0 for i in range(0,11,10)]
    run_scan(systems,points,scan_range,corr_torsion_scan,"Torsion bend",tester)

    ###Improper torsion scan of ammonia H-NH-H angle
    ammonia=corr_imp_scan[0].as_universe()
    systems=[psr.System(ammonia,False) for i in range(2)]
    for i in range(2):
        systems[0].insert(ammonia[i])
        systems[1].insert(ammonia[i+2])
    points=[ammonia[i] for i in range(4)]
    scan_range=[i*math.pi/180.0 for i in range(0,11,10)]
    run_scan(systems,points,scan_range,corr_imp_scan,"Imp bend",tester)

    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
