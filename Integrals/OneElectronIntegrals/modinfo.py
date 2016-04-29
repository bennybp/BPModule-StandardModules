from pulsar.datastore import OptionType

minfo = {

  "Overlap" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "OneElectronIntegrals.so",
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
    "modpath"     : "OneElectronIntegrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO kinetic energy integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                    }
  },

}


