from pulsar.datastore import OptionType
from pulsar.datastore.OptionValidators import *

MethodOption=(OptionType.String,None,True,None,
                    'The key of the method that should be used to compute the '\
                    'energy derivative.')
FraggerOption=(OptionType.String,"PSR_BOND_FRAG",False,None,
               'The key used to fragment the system')

DerivOption=(OptionType.Int,255,False,None,
                    'What order analytic derivatives are available (actually '\
                    'have arbitrary order available)')

BasisOption=(OptionType.String,"Primary",False,None,
                    'What basis set tag should be used')
MIMOption=(OptionType.String,"PSR_MIM",False,None,
                    "A way for changing which MIM module is called.")

minfo = {
  "CP":{
    "type"        : "python_module",
    "base"        : "EnergyMethod",
    "version"     : "0.1a",
    "description" : "Performs a standard CP correctio on a system",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {"METHOD" : MethodOption,
                     "GHOSTER_KEY": (OptionType.String,"PSR_CP_FRAG",False,None,
                     'The key used to put ghosts on the system'),
                     "MBE_KEY":(OptionType.String,"PSR_MBE",False,None,
                     'The key used for MBEs on the system.  For CP,'\
                     ' the MBE should be truncated at order 1')
                    }
  },
  "MBE" :
  {
    "type"        : "c_module",
    "base"        : "EnergyMethod",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Performs a many-body expansion on a system",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {"METHOD" : MethodOption,
                     "FRAGMENTIZER": FraggerOption
                    }
  },
  "MIM" :
  {
    "type"        : "c_module",
    "base"        : "EnergyMethod",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Performs a Molecules in Molecules (MIM) computation",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "WEIGHTS":(OptionType.ListFloat,None,True,None,
                    'WEIGHTS[i] is the i-th weight of the i-th system using '\
                    'the i-th method.'),
                    "METHODS":(OptionType.ListString,None,True,None,
                    "METHODS[i] is the i-th method's key, you may provide only"\
                    ' one key if it is the systems that are changing'),
                    "FRAGMENTIZER":(OptionType.String,"BP_NULL_FRAG",False,None,
                    'The key used to fragment the system'),
                    }
  },
  "SCF" :
  {
    "type"         : "python_module",
    "base"        : "EnergyMethod",
    "version"     : "0.1a",
    "description" : "Calls Psi4 via a system call and then runs an SCF",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""], 
    "options"     : {
                    }
  },
  "MP2" :
  {
    "type"        : "python_module",
    "base"        : "EnergyMethod",
    "version"     : "0.1a",
    "description" : "Calls Psi4 via a system call and then runs MP2",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""], 
    "options"     : {
                    }
  },
  "CCSD(T)" :
  {
    "type"        : "python_module",
    "base"        : "EnergyMethod",
    "version"     : "0.1a",
    "description" : "Calls Psi4 via a system call and then runs CCSD(T)",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""], 
    "options"     : {
                    }
  },
  "HelgakerCBS" :
    {
    "type"        : "python_module",
    "base"        : "EnergyMethod",
    "version"     : "0.1a",
    "description" : "Performs a Complete Basis Set Extrapolation using the two"\
                    "-point Helgaker formula",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "BASIS_CARDINAL_NUMS":(OptionType.ListInt,None,True,None,
                    "The cardinal numbers of the two basis sets."),
                    "MIM_KEY":MIMOption
                    }
  },
  "FellerCBS" :
    {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Performs a Complete Basis Set Extrapolation using the "\
                    "three-point Feller formula",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "MAX_DERIV":(OptionType.Int,2,False,None,"The maximum "\
                    "analytic derivative available.  At the moment this is 2"),
                    "BASIS_CARDINAL_NUMS":(OptionType.ListInt,None,True,None,
                    "The cardinal numbers of the two basis sets."),
                    "MIM_KEY":MIMOption
                    }
  },
  "FPA" :
    {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Performs a Focal Point Analysis",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "MAX_DERIV":DerivOption,
                    "LARGE_MP2_KEY":(OptionType.String,None,True,None,
                       "The key for the large MP2 module"),
                    "SMALL_MP2_KEY":(OptionType.String,"BP_MP2",False,None,
                       "The key for the small MP2 module"),
                    "CCSD(T)_KEY":(OptionType.String,"BP_CCSD(T)",False,None,
                       "The key for the CCSD(T) module"),
                    "MIM_KEY":MIMOption
                    }
  },
  "HFIterate" :
  {
    "type"        : "c_module",
    "base"        : "SCFIterator",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Quick HF test calculation",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [""],
    "options"     : {
                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
                            "Key of the ao overlap module to use"),
                        "BASIS_SET"       :  (OptionType.String, "Primary", False, None,
                            'Tag representing the basis set in the system'),
                        "KEY_AO_CACHER": (OptionType.String, None, True, None,
                            "Key of the one-electron integral cacher"),
                    }
  },
  "BasicFockBuild" :
  {
    "type"        : "c_module",
    "base"        : "FockBuilder",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Quick HF test calculation",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [""],
    "options"     : {
                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
                            "Key of the ao overlap module to use"),
                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
                            "Key of the core builder module to use"),
                        "KEY_AO_CACHER": (OptionType.String, None, True, None,
                            "Key of the one-electron integral cacher"),
                        "KEY_AO_ERI": (OptionType.String, None, True, None,
                            "Key of the ERI module to use"),
                    }
  },
  "Damping" :
  {
    "type"        : "c_module",
    "base"        : "EnergyMethod",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Quick HF test calculation",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [""],
    "options"     : {
                        "KEY_INITIAL_GUESS": (OptionType.String, None, False, None,
                            "Key for the initial guess module"),
                        "KEY_SCF_ITERATOR": (OptionType.String, None, True, None,
                            "Key of the iterator module to use"),
                        "KEY_FOCK_BUILDER": (OptionType.String, None, True, None,
                            "Key of the fock builder module to use"),
                        "MAX_ITER": (OptionType.Int, 40, False, None,
                            "Key of the ao electron repulsion integral module to use"),
                        "DENS_TOLERANCE": (OptionType.Float, 1e-8, False, None,
                            "Maximum value for the change in density"),
                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
                            "Key of the core builder module to use"),
                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
                            "Key of the nuclear repulsion module to use"),
                        "KEY_AO_CACHER": (OptionType.String, None, True, None,
                            "Key of the one-electron integral cacher"),
                        "DAMPING_FACTOR": (OptionType.Float, 0.0, False, RangeCheck(0.0, 1.0, True, False),
                            "Amount of old fock matrix to use in constructing new fock matrix"),
                    }
  },
  "DIIS" :
  {
    "type"        : "c_module",
    "base"        : "EnergyMethod",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Quick HF test calculation",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [""],
    "options"     : {
                        "KEY_INITIAL_GUESS": (OptionType.String, None, False, None,
                            "Key for the initial guess module"),
                        "KEY_SCF_ITERATOR": (OptionType.String, None, True, None,
                            "Key of the iterator module to use"),
                        "KEY_FOCK_BUILDER": (OptionType.String, None, True, None,
                            "Key of the fock builder module to use"),
                        "MAX_ITER": (OptionType.Int, 40, False, None,
                            "Key of the ao electron repulsion integral module to use"),
                        "DENS_TOLERANCE": (OptionType.Float, 1e-8, False, None,
                            "Maximum value for the change in density"),
                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
                            "Key of the core builder module to use"),
                        "KEY_AO_CACHER": (OptionType.String, None, True, None,
                            "Key of the one-electron integral cacher"),
                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
                            "Key of the nuclear repulsion module to use"),
                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
                            "Key of the ao overlap module to use"),
                    }
  },
  "CoreGuess" :
  {
    "type"        : "c_module",
    "base"        : "EnergyMethod",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Initial guess for Hartree Fock via core guess",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [""],
    "options"     : {
                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
                            "Key of the nuclear repulsion module to use"),
                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
                            "Key of the ao overlap module to use"),
                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
                            "Key of the core builder module to use"),
                        "KEY_AO_CACHER": (OptionType.String, None, True, None,
                            "Key of the one-electron integral cacher"),
                    }
  },
}


