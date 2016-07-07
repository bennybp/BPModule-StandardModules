from pulsar.datastore import OptionType
from pulsar.datastore.OptionValidators import *

minfo = {

  "OSOverlap" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO overlap integrals over gaussian basis functions via Obara-Saika",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "OSDipole" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO overlap integrals over gaussian basis functions via Obara-Saika",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "OSKineticEnergy" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO kinetic energy integrals over gaussian basis functions via Obara-Saika",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "OSOneElectronPotential" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO electron-nuclear attraction integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "GRID":   ( OptionType.String,  None, True, None,  "Grid of point charges to calculate the potential with" )
                    }
  },

  "OneElectronIntegralSum" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Building of sum of one-electron integrals",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "KEY_AO_TERMS":   ( OptionType.ListString,  None, True, None,  "Keys to the one-electron integral modules to sum")
                    }
  },
  "OneElectronProperty" :
  {
    "type"        : "c_module",
    "base"        : "PropertyCalculator",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "A general one-electron property calculator",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "KEY_ONEEL_MOD":   ( OptionType.String,  None, True, None,  "Key of which one electron integral to use"),
                    }
  },
  "OneElectron_Eigen" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronMatrix",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Caching of one-electron integrals in an Eigen3 matrix",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "CACHE_RESULTS":   ( OptionType.Bool,  True, False, None,  "Cache the results between instantiations"),
                    }
  },

  "ReferenceERI" :
  {
    "type"        : "c_module",
    "base"        : "TwoElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation ERI via a slow but accurate method",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "NuclearRepulsion" :
  {
    "type"        : "c_module",
    "base"        : "SystemIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of nuclear-nuclear repulsion",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },
  "NuclearDipole" :
  {
    "type"        : "c_module",
    "base"        : "SystemIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of nuclear dipole moment",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },
}


