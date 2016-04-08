from bpmodule.datastore import OptionType

DistThreshOption=(OptionType.DictIntFloat,{},False,
                    None,'A int->float dictionary where the int, call it n, is'\
                    ' the n-mer and the float is the maximum distance '\
                    'that the n monomers can be apart, e.g. {2:3.14} means '\
                    'that dimers whose centers of mass are more than 3.14 '\
                    'Angstroms apart are excluded.')   

TruncOrderOption=(OptionType.Int,1,False,None,
                 'Unions of up to how many fragments should be created?'\
                 'In other words, do you want dimers, trimers,... to be made'\
                 'as well?')

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
                    "TRUNCATION_ORDER": TruncOrderOption,
                    "DISTANCE_THRESHOLDS": DistThreshOption
                    }
  },
  "UserDefined" :
  {
    "type"        : "c_module",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes all atoms within N bonds a fragment",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                     "FRAGMENT_NAMES":(OptionType.ListString,None,True,None,
                     'Names of the fragments'),
                     "ATOMS_PER_FRAG":(OptionType.ListInt,None,True,None,
                     'Number of atoms in each fragment'),
                     "FRAGMENTS":(OptionType.ListInt,None,True,None,
                     'The atoms in each fragment'),
                     "TRUNCATION_ORDER": TruncOrderOption,
                     "DISTANCE_THRESHOLDS":DistThreshOption 
                    }
  },
  "NullFragmenter" :
  {
    "type"        : "c_module",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Returns the input as a single fragment.",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                      "TRUNCATION_ORDER":(OptionType.Int,0,False,None,
                      "Should not be set for a NullFragmenter"),
                      "DISTANCE_THRESHOLDS": DistThreshOption
                    }
  },
}


