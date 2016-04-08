import sys
import argparse
import traceback
import bpmodule as bp

bp.Init(sys.argv, out = "stdout", color = True, debug = True)

def LoadDefaultModules(mm):
   mm.LoadModule("Methods","CP","BP_CP")
   mm.LoadModule("Methods","VMFC","BP_VMFC")
   mm.LoadModule("Methods","MIM","BP_MIM")
   mm.LoadModule("Methods","MBE","BP_MBE")
   mm.LoadModule("Methods","SCF","BP_SCF")
   mm.LoadModule("Methods","HelgakerCBS","BP_HELGAKER_CBS")
   mm.LoadModule("SystemFragmenters","UserDefined","BP_UD_FRAG")
   mm.LoadModule("SystemFragmenters","NullFragmenter","BP_NULL_FRAG")
   mm.LoadModule("SystemFragmenters","Atomizer","BP_ATOM_FRAG")
   mm.LoadModule("SystemFragmenters","Bondizer","BP_BOND_FRAG")
