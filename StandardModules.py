import sys
import argparse
import traceback
import pulsar as psr

psr.initialize(sys.argv, out = "stdout", color = True, debug = True)

def LoadDefaultModules(mm):
   #mm.load_module("Methods","MIM","PSR_MIM")
   mm.load_module("Methods","MBE","PSR_MBE")
   mm.load_module("Methods","CP","PSR_CP")
   mm.load_module("Methods","FPA","PSR_FPA")
   mm.load_module("Methods","HelgakerCBS","PSR_HELGAKER_CBS")
   mm.load_module("Methods","FellerCBS","PSR_FELLER_CBS")
   mm.load_module("SystemFragmenters","UserDefined","PSR_UD_FRAG")
   mm.load_module("SystemFragmenters","NullFragmenter","PSR_NULL_FRAG")
   mm.load_module("SystemFragmenters","Atomizer","PSR_ATOM_FRAG")
   mm.load_module("SystemFragmenters","Bondizer","PSR_BOND_FRAG")
   mm.load_module("SystemFragmenters","Ghoster","PSR_GHOST_FRAG")
   mm.load_module("SystemFragmenters","CPGhoster","PSR_CP_FRAG")
   mm.load_module("SystemFragmenters","VMFCGhoster","PSR_VMFC_FRAG")
   mm.load_module("SystemFragmenters","CrystalFragger","PSR_CRYS_FRAG")
