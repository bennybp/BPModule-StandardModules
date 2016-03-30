from bpmodule.datastore import OptionType

minfo = {

  "MBE" :
  {
    "type"        : "c_module",
    "modpath"     : "Methods.so",
    "version"     : "0.1a",
    "description" : "Performs a many-body expansion on a system",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "METHOD":(OptionType.String,None,True,None,
                    'The key of the method that should be used to compute the '\
                    'energy derivative.'),
                    "FRAGMENTIZER":(OptionType.String,"FRAG",False,None,
                    'The key used to fragment the system'),
                    }
  },
  "MIM" :
  {
    "type"        : "c_module",
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
                    "FRAGMENTIZER":(OptionType.String,"FRAG",False,None,
                    'The key used to fragment the system'),
                    }
  },
  "SCF" :
  {
    "type"         : "python_module",
    "version"     : "0.1a",
    "description" : "Calls Psi4 via a system call and then runs an SCF",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""], 
    "options"     : {}
  },

}


