from pulsar.datastore import OptionType

TruncOrderOption=(OptionType.Int,1,False,None,
        'Unions of up to how many fragments should be created?'\
        'In other words, do you want dimers, trimers,... to be made as well?')
SubFragger=(OptionType.String,"PSR_BOND_FRAG",False,None,
        'The SystemFragmenter used to chop the supercell into fragments.')
GhosterKey=(OptionType.String,"PSR_GHOST_FRAG",False,None,
        'The generic ghosted frament maker.')
    
minfo = {

  "Atomizer" :
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
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
    "base"        : "SystemFragmenter",
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
                    }
  },
  "CrystalFragger" :
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Fragments a periodic system",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                      "SYSTEM_FRAGMENTER_KEY":SubFragger,
                      "LATTICE_SIZE":(OptionType.ListInt,[3,3,3],False,None,
                      'How large of a supercell should we fragment?'),
                    },
  },
  #TODO: Can I get rid of this middle fragmenter, if not why?
  "Ghoster":
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes fragments with ghost atoms",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                      "SYSTEM_FRAGMENTER_KEY":SubFragger,
                      "GHOST_TRUNCATION_ORDERS":(OptionType.DictIntInt,{},False,
                      None,'A map describing the maximum number of ghost'\
                      ' fragments, g, for a real-mer.  The map is of the form '
                      ' {{r:g}}.  Any r not specified is assumed to be 0.')
                    }
  },
  "CPGhoster":
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes fragments with ghost atoms consistent with CP",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                     "GHOSTER_KEY":GhosterKey,
                    }
  }, 
  "VMFCGhoster":
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes fragments with ghost atoms consistent with VMFC",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                     "GHOSTER_KEY":GhosterKey
                    }
  },
  "UserDefined" :
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Makes whatever you asked for into fragments",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                     "FRAGMENT_NAMES":(OptionType.ListString,None,True,None,
                     'Names of the fragments'),
                     "ATOMS_PER_FRAG":(OptionType.ListInt,None,True,None,
                     'Number of atoms in each fragment'),
                     "FRAGMENTS":(OptionType.ListInt,None,True,None,
                     'The atoms in each fragment'),
                    }
  },
  "NullFragmenter" :
  {
    "type"        : "c_module",
    "base"        : "SystemFragmenter",
    "modpath"     : "SystemFragmenters.so",
    "version"     : "0.1a",
    "description" : "Returns the input as a single fragment.",
    "authors"     : ["Ryan Richard"],
    "refs"        : [""],
    "options"     : {
                      "TRUNCATION_ORDER":(OptionType.Int,0,False,None,
                      "Should not be set for a NullFragmenter"),
                    }
  },
}


