from pulsar.datastore import OptionType
from pulsar.datastore.OptionValidators import *

minfo = {

  "Overlap" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO overlap integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "Dipole" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO overlap integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "KineticEnergy" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO kinetic energy integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

  "OneElectronPotential" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO nuclear repulsion integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "GRID":   ( OptionType.String,  None, True, None,  "Grid of point charges to calculate the potential with" )
                    }
  },

  "CoreBuild" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Building of the core hamiltonian",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "KEY_AO_CORE_TERMS":   ( OptionType.ListString,  None, True, None,  "One-electron integrals to use in the core")
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
  "EigenCacher" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronCacher",
    "modpath"     : "Integrals.so",
    "version"     : "0.1a",
    "description" : "Caching of one-electron integrals in an Eigen3 matrix",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
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


