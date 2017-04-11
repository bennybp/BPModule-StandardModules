from pulsar import OptionType

#Stuff common to many modules
c_mod,sysfrag,modpath="c_module","SystemFragmenter","pulsar_modules.so"
version="0.1a"
no_options={}
no_ref,Ben,Ryan=[""],["Ben Pritchard"],["Ryan Richard"]

sys_frags=["Atomizer","Bondizer","CrystalFragger","CPGhoster","VMFCGhoster",\
           "NMerizer","MBE"]
minfo={i:{
    "type":c_mod,
    "base":sysfrag,
    "modpath":modpath,
    "version":version,
    "authors":Ryan,
    "refs":no_ref,
    "options":no_options} 
for i in sys_frags}

minfo["Atomizer"]["description"]="Atomizes a System"
minfo["Atomizer"]["authors"]=Ben

minfo["Bondizer"]["description"]="Makes all atoms within N bonds a fragment"
minfo["Bondizer"]["options"]={
                    "MAX_NBONDS":(OptionType.Int,2147483647,False,None,
                    'The max number of bonds between atoms, by default all'\
                    ' bonded atoms are in a fragment.'),
                    }

minfo["CrystalFragger"]["description"]="Fragments a periodic system"
#sf["crystal"]["options"]={
#                      "SYSTEM_FRAGMENTER_KEY":SubFragger,
#                      "LATTICE_SIZE":(OptionType.ListInt,[3,3,3],False,None,
#                      'How large of a supercell should we fragment?'),
#                    }
minfo["CPGhoster"]["description"]="Makes fragments with ghost atoms consistent with CP"
minfo["CPGhoster"]["options"]={"SYSTEM_FRAGMENTER_KEY":(OptionType.String,None,True,None,
                            "Fragmenter to call for original fragments")}
minfo["VMFCGhoster"]["description"]="Makes fragments with ghost atoms consistent with VMFC"
minfo["VMFCGhoster"]["options"]={"SYSTEM_FRAGMENTER_KEY":(OptionType.String,None,True,
                              None,"Fragmenter to call to get frags")}
minfo["NMerizer"]["description"]="Makes fragments from union of subfragments"
minfo["NMerizer"]["options"]={
                    "SYSTEM_FRAGMENTER_KEY":(OptionType.String,None,True,None,
                    "The key to generate the first set of fragments"),
                    "TRUNCATION_ORDER":(OptionType.Int,1,False,None,
                    "The maximum number of fragments involved in a union"),
                    "DISTANCE_THRESHOLDS":(OptionType.DictIntFloat,{},False,
                     None,"Maximum distance per truncation order to use") 
                    }
minfo["MBE"]["base"]="EnergyMethod"
minfo["MBE"]["description"]="Runs a many-body expansion"
minfo["MBE"]["options"]={
                "SYSTEM_FRAGMENTER_KEY":(OptionType.String,None,True,None,
                "Fragmenter to call for original fragments"),
                "METHOD_KEY":(OptionType.String,None,True,None,
                "EnergyMethod to call"),
                    }

#  "CP":{
#    "type"        : "python_module",
#    "base"        : "EnergyMethod",
#    "version"     : "0.1a",
#    "description" : "Performs a standard CP correctio on a system",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {"METHOD" : MethodOption,
#                     "GHOSTER_KEY": (OptionType.String,"PSR_CP_FRAG",False,None,
#                     'The key used to put ghosts on the system'),
#                     "MBE_KEY":(OptionType.String,"PSR_MBE",False,None,
#                     'The key used for MBEs on the system.  For CP,'\
#                     ' the MBE should be truncated at order 1')
#                    }
#  },
#  "MBE" :
#  {
#    "type"        : "c_module",
#    "base"        : "EnergyMethod",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Performs a many-body expansion on a system",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {"METHOD" : MethodOption,
#                     "FRAGMENTIZER": FraggerOption
#                    }
#  },
#  "CorrelationEnergy" :
#    {
#    "type"        : "python_module",
#    "base"        : "EnergyMethod",
#    "version"     : "0.1a",
#    "description" : "Computes correlation energy by subtacting refernce",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {
#                    "CORRELATED_KEY":(OptionType.String,None,True,None,
#                    "The key for the method that generates the correlated egy"),
#                    "REFERENCE_KEY":(OptionType.String,None,True,None,
#                    "The key for the method that generates the reference egy")
#                    }
#  },
#  "HelgakerCBS" :
#    {
#    "type"        : "python_module",
#    "base"        : "EnergyMethod",
#    "version"     : "0.1a",
#    "description" : "Performs a Complete Basis Set Extrapolation using the two"\
#                    "-point Helgaker formula",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {
#                    "BASIS_CARDINAL_NUMS":(OptionType.ListInt,None,True,None,
#                    "The cardinal numbers of the two basis sets."),
#                    "METHODS":(OptionType.ListString,None,True,None,
#                    "The respective methods for the cardinal numbers"
#                    )
#                    }
#  },
#  "FellerCBS" :
#    {
#    "type"        : "python_module",
#    "version"     : "0.1a",
#    "description" : "Performs a Complete Basis Set Extrapolation using the "\
#                    "three-point Feller formula",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {
#                    "MAX_DERIV":(OptionType.Int,2,False,None,"The maximum "\
#                    "analytic derivative available.  At the moment this is 2"),
#                    "BASIS_CARDINAL_NUMS":(OptionType.ListInt,None,True,None,
#                    "The cardinal numbers of the two basis sets."),
#                    "METHODS":(OptionType.ListString,None,True,None,
#                    "Keys of methods for the respective cardinal numbers"
#                    )
#                    }
#  },
#  "FPA" :
#    {
#    "type"        : "python_module",
#    "version"     : "0.1a",
#    "description" : "Performs a Focal Point Analysis",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {
#                    "MAX_DERIV":DerivOption,
#                    "LARGE_MP2_KEY":(OptionType.String,None,True,None,
#                       "The key for the large MP2 module"),
#                    "SMALL_MP2_KEY":(OptionType.String,"BP_MP2",False,None,
#                       "The key for the small MP2 module"),
#                    "CCSD(T)_KEY":(OptionType.String,"BP_CCSD(T)",False,None,
#                       "The key for the CCSD(T) module"),
#                    "MIM_KEY":MIMOption
#                    }
#  },
#  "MyCrzyCompMeth" :
#    {
#    "type"        : "python_module",
#    "version"     : "0.1a",
#    "description" : "Proof of concept crzy nested method",
#    "authors"     : ["Ryan Richard"],
#    "refs"        : [""],
#    "options"     : {
#                    "MAX_DERIV":(OptionType.Int,2,False,None,"The maximum "\
#                    "analytic derivative available.  At the moment this is 2"),
#                    }
#  },
#  "HFIterate" :
#  {
#    "type"        : "c_module",
#    "base"        : "SCFIterator",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Quick HF test calculation",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [""],
#    "options"     : {
#                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
#                            "Key of the ao overlap module to use"),
#                        "BASIS_SET"       :  (OptionType.String, "Primary", False, None,
#                            'Tag representing the basis set in the system'),
#                        "KEY_ONEEL_MAT": (OptionType.String, None, True, None,
#                            "Key of the one-electron integral cacher"),
#                    }
#  },
#  "BasicFockBuild" :
#  {
#    "type"        : "c_module",
#    "base"        : "FockBuilder",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Quick HF test calculation",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [""],
#    "options"     : {
#                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
#                            "Key of the ao overlap module to use"),
#                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
#                            "Key of the core builder module to use"),
#                        "KEY_ONEEL_MAT": (OptionType.String, None, True, None,
#                            "Key of the one-electron integral cacher"),
#                        "KEY_AO_ERI": (OptionType.String, None, True, None,
#                            "Key of the ERI module to use"),
#                    }
#  },
#  "Damping" :
#  {
#    "type"        : "c_module",
#    "base"        : "EnergyMethod",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Quick HF test calculation",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [""],
#    "options"     : {
#                        "KEY_INITIAL_GUESS": (OptionType.String, None, False, None,
#                            "Key for the initial guess module"),
#                        "KEY_SCF_ITERATOR": (OptionType.String, None, True, None,
#                            "Key of the iterator module to use"),
#                        "KEY_FOCK_BUILDER": (OptionType.String, None, True, None,
#                            "Key of the fock builder module to use"),
#                        "MAX_ITER": (OptionType.Int, 40, False, None,
#                            "Key of the ao electron repulsion integral module to use"),
#                        "DENS_TOLERANCE": (OptionType.Float, 1e-8, False, None,
#                            "Maximum value for the change in density"),
#                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
#                            "Key of the core builder module to use"),
#                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
#                            "Key of the nuclear repulsion module to use"),
#                        "KEY_ONEEL_MAT": (OptionType.String, None, True, None,
#                            "Key of the one-electron integral cacher"),
#                        "DAMPING_FACTOR": (OptionType.Float, 0.0, False, RangeCheck(0.0, 1.0, True, False),
#                            "Amount of old fock matrix to use in constructing new fock matrix"),
#                    }
#  },
#  "DIIS" :
#  {
#    "type"        : "c_module",
#    "base"        : "EnergyMethod",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Quick HF test calculation",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [""],
#    "options"     : {
#                        "KEY_INITIAL_GUESS": (OptionType.String, None, False, None,
#                            "Key for the initial guess module"),
#                        "KEY_SCF_ITERATOR": (OptionType.String, None, True, None,
#                            "Key of the iterator module to use"),
#                        "KEY_FOCK_BUILDER": (OptionType.String, None, True, None,
#                            "Key of the fock builder module to use"),
#                        "MAX_ITER": (OptionType.Int, 40, False, None,
#                            "Key of the ao electron repulsion integral module to use"),
#                        "DENS_TOLERANCE": (OptionType.Float, 1e-8, False, None,
#                            "Maximum value for the change in density"),
#                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
#                            "Key of the core builder module to use"),
#                        "KEY_ONEEL_MAT": (OptionType.String, None, True, None,
#                            "Key of the one-electron integral cacher"),
#                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
#                            "Key of the nuclear repulsion module to use"),
#                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
#                            "Key of the ao overlap module to use"),
#                    }
#  },
#  "CoreGuess" :
#  {
#    "type"        : "c_module",
#    "base"        : "EnergyMethod",
#    "modpath"     : "Methods.so",
#    "version"     : "0.1a",
#    "description" : "Initial guess for Hartree Fock via core guess",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [""],
#    "options"     : {
#                        "KEY_NUC_REPULSION": (OptionType.String, None, True, None,
#                            "Key of the nuclear repulsion module to use"),
#                        "KEY_AO_OVERLAP": (OptionType.String, None, True, None,
#                            "Key of the ao overlap module to use"),
#                        "KEY_AO_COREBUILD": (OptionType.String, None, True, None,
#                            "Key of the core builder module to use"),
#                        "KEY_ONEEL_MAT": (OptionType.String, None, True, None,
#                            "Key of the one-electron integral matrix generator"),
#                    }
#  },
#  "OSOverlap" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of AO overlap integrals over gaussian basis functions via Obara-Saika",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#
#  "OSDipole" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of AO overlap integrals over gaussian basis functions via Obara-Saika",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#
#  "OSKineticEnergy" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of AO kinetic energy integrals over gaussian basis functions via Obara-Saika",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#
#  "OSOneElectronPotential" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of AO electron-nuclear attraction integrals over gaussian basis functions",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                        "GRID":   ( OptionType.String,  None, True, None,  "Grid of point charges to calculate the potential with" )
#                    }
#  },
#
#  "OneElectronIntegralSum" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Building of sum of one-electron integrals",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                        "KEY_AO_TERMS":   ( OptionType.ListString,  None, True, None,  "Keys to the one-electron integral modules to sum")
#                    }
#  },
#  "OneElectronProperty" :
#  {
#    "type"        : "c_module",
#    "base"        : "PropertyCalculator",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "A general one-electron property calculator",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                        "KEY_ONEEL_MOD":   ( OptionType.String,  None, True, None,  "Key of which one electron integral to use"),
#                    }
#  },
#  "OneElectron_Eigen" :
#  {
#    "type"        : "c_module",
#    "base"        : "OneElectronMatrix",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Caching of one-electron integrals in an Eigen3 matrix",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                        "CACHE_RESULTS":   ( OptionType.Bool,  True, False, None,  "Cache the results between instantiations"),
#                    }
#  },
#
#  "ReferenceERI" :
#  {
#    "type"        : "c_module",
#    "base"        : "TwoElectronIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation ERI via a slow but accurate method",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#
#  "NuclearRepulsion" :
#  {
#    "type"        : "c_module",
#    "base"        : "SystemIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of nuclear-nuclear repulsion",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#  "NuclearDipole" :
#  {
#    "type"        : "c_module",
#    "base"        : "SystemIntegral",
#    "modpath"     : "Integrals.so",
#    "version"     : "0.1a",
#    "description" : "Calculation of nuclear dipole moment",
#    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
#    "refs"        : [],
#    "options"     : {
#                    }
#  },
#}
#
#
#
