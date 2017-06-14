from .pulsar_modules import *
from .modinfo import *
from .creator import *

def initialize(mm):
    mm.load_module("pulsar_modules","Atomizer","PSR_ATOM_FRAG")
    mm.load_module("pulsar_modules","NMerizer","PSR_NMER_FRAG")
    mm.load_module("pulsar_modules","Bondizer","PSR_BOND_FRAG")
    mm.load_module("pulsar_modules","CPGhoster","PSR_CP_FRAG")
    mm.load_module("pulsar_modules","VMFCGhoster","PSR_VMFC_FRAG")
    mm.load_module("pulsar_modules","MBE","PSR_MBE")
