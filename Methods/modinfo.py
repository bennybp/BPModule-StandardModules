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
                    "TRUNCATION_ORDER":(OptionType.Int,3,False,None,
                    "The order at which the MBE should be truncated"),
                    "DISTANCE_THRESHOLDS":(OptionType.DictIntFloat,{},False,
                    None,'A int->float dictionary where the int, call it n, is'\
                    ' the MBE order and the float is the maximum distance '\
                    'that the n monomers can be apart, e.g. {2:3.14} means '\
                    'that dimers whose centers of mass are more than 3.14 '\
                    'Angstroms apart are excluded.')
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


