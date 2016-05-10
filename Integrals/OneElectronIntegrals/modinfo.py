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

  "OneElectronPotential" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "OneElectronIntegrals.so",
    "version"     : "0.1a",
    "description" : "Calculation of AO nuclear repulsion integrals over gaussian basis functions",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "grid":   ( OptionType.String,  None, True, None,  "Grid of point charges to calculate the potential with" )
                    }
  },

  "CoreBuild" :
  {
    "type"        : "c_module",
    "base"        : "OneElectronIntegral",
    "modpath"     : "OneElectronIntegrals.so",
    "version"     : "0.1a",
    "description" : "Building of the core hamiltonian",
    "authors"     : ["Benjamin Pritchard <ben@bennyp.org>"],
    "refs"        : [],
    "options"     : {
                        "KEY_AO_KINETIC":   ( OptionType.String,  None, True, None,  "Key of which kinetic energy module to use"),
                        "KEY_AO_NUCATT":   ( OptionType.String,  None, True, None,  "Key of which electron-nuclear attraction module to use"),
                        "KEY_AO_ADDITIONAL":   ( OptionType.ListString,  None, False, None,  "Additional one-electron integrals to incorporate ")
                    }
  },

}


