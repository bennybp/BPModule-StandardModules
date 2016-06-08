import sys
import argparse
import traceback
import pulsar as psr

psr.Init(sys.argv, out = "stdout", color = True, debug = True)

def LoadDefaultModules(mm):
   #mm.LoadModule("Methods","MIM","PSR_MIM")
   mm.LoadModule("Methods","MBE","PSR_MBE")
   mm.LoadModule("Methods","SCF","PSR_SCF")
   mm.LoadModule("Methods","MP2","PSR_MP2")
   mm.LoadModule("Methods","CCSD(T)","PSR_CCSD(T)")
   mm.LoadModule("Methods","CP","PSR_CP")
   mm.LoadModule("Methods","FPA","PSR_FPA")
   mm.LoadModule("Methods","HelgakerCBS","PSR_HELGAKER_CBS")
   mm.LoadModule("Methods","FellerCBS","PSR_FELLER_CBS")
   mm.LoadModule("SystemFragmenters","UserDefined","PSR_UD_FRAG")
   mm.LoadModule("SystemFragmenters","NullFragmenter","PSR_NULL_FRAG")
   mm.LoadModule("SystemFragmenters","Atomizer","PSR_ATOM_FRAG")
   mm.LoadModule("SystemFragmenters","Bondizer","PSR_BOND_FRAG")
   mm.LoadModule("SystemFragmenters","Ghoster","PSR_GHOST_FRAG")
   mm.LoadModule("SystemFragmenters","CPGhoster","PSR_CP_FRAG")
   mm.LoadModule("SystemFragmenters","VMFCGhoster","PSR_VMFC_FRAG")
   mm.LoadModule("SystemFragmenters","CrystalFragger","PSR_CRYS_FRAG")
