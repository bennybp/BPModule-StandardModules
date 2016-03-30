from bpmodule.datastore import OptionType

minfo = {

  "Atomizer" :
  {
    "type"        : "c_module",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Atomizes a System",
    "authors"     : ["Ben Pritchard"],
    "refs"        : [""],
    "options"     : {
                        # No options
                    }
  },
  "Bondizer" :
  {
    "type"        : "c_module",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes all atoms within N bonds a fragment",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                    "NBONDS":(OptionType.Int,2147483647,False,None,
                    'The max number of bonds between atoms for them to be '\
                    'considered in the same fragment.  By default all '\
                    'covalently bonded atoms are in the same fragment '\
                    '(assuming they are within 2,147,483,647 bonds of each '\
                    'other).'),
                    "TRUNCATION_ORDER":(OptionType.Int,1,False,None,
                    'Unions of up to how many fragments should be created?'\
                    'In other words, do you want dimers, trimers,... to be made'\
                    'as well?'),
                    "DISTANCE_THRESHOLDS":(OptionType.DictIntFloat,{},False,
                    None,'A int->float dictionary where the int, call it n, is'\
                    ' the n-mer and the float is the maximum distance '\
                    'that the n monomers can be apart, e.g. {2:3.14} means '\
                    'that dimers whose centers of mass are more than 3.14 '\
                    'Angstroms apart are excluded.')   
                    }
  },

}


