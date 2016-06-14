from pulsar.datastore import OptionType

minfo = {

  "NuclearRepulsion" :
  {
    "type"        : "c_module",
    "base"        : "SystemIntegral",
    "modpath"     : "SystemIntegrals.so",
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
    "modpath"     : "SystemIntegrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of nuclear dipole moment",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

}


