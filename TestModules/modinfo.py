from bpmodule.datastore import OptionType

minfo = {

  "TestModule1" :
  {
    "type"        : "c_module",
    "modpath"     : "TestModules.so",
    "version"     : "0.1a",
    "description" : "Some test module",
    "authors"     : ["me", "myself", "I"],
    "refs"        : ["some paper", "some other paper"],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                        "bool_opt_def":    (  OptionType.Bool,  False,     False,   None,  "Some help string"),
                        "bool_opt":        (  OptionType.Bool,  None,      False,   None,  "Some help string"),

                        "int_opt_def":     (  OptionType.Int,    1,        False,   None,  "Some help string"),
                        "int_opt":         (  OptionType.Int,    None,     False,   None,  "Some help string"),

                        "double_opt_def":  ( OptionType.Float,   6.02e23,  False,   None,  "Some help string"),
                        "double_opt":      ( OptionType.Float,   None,     False,   None,  "Some help string"),

                        "str_opt_def":     (  OptionType.String,    "Hi",     False,   None,  "Some help string"),
                        "str_opt":         (  OptionType.String,    None,     False,   None,  "Some help string"),
                    }
  },


 "TestPyModule1" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Some test module",
    "authors"     : ["me", "myself", "I"],
    "refs"        : ["some paper", "some other paper"],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                        "bool_opt_def":    (  OptionType.Bool,  False,     False,   None,  "Some help string"),
                        "bool_opt":        (  OptionType.Bool,  None,      False,   None,  "Some help string"),

                        "int_opt_def":     (  OptionType.Int,    1,        False,   None,  "Some help string"),
                        "int_opt":         (  OptionType.Int,    None,     False,   None,  "Some help string"),

                        "double_opt_def":  ( OptionType.Float,   6.02e23,  False,   None,  "Some help string"),
                        "double_opt":      ( OptionType.Float,   None,     False,   None,  "Some help string"),

                        "str_opt_def":     (  OptionType.String,    "Hi",     False,   None,  "Some help string"),
                        "str_opt":         (  OptionType.String,    None,     False,   None,  "Some help string"),
                    }
  },


  ##########################################
  # Testing of linking to external libraries
  ##########################################
  "TestExtLib" :
  {
    "type"        : "c_module",
    "modpath"     : "TestModules.so",
    "version"     : "0.1a",
    "description" : "Some test module",
    "authors"     : ["me", "myself", "I"],
    "refs"        : ["some paper", "some other paper"],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                    }
  },



  ##########################################
  # Testing of options parsing. These modules
  # aren't meant to be loaded
  ##########################################
  "TestOptions_int" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of integers",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                        "int_opt_def":     (  OptionType.Int,   5,         False,   None,  "Some help string"),
                        "int_req":         (  OptionType.Int,   None,      True,    None,  "Some help string"),
                        "int_opt":         (  OptionType.Int,   None,      False,   None,  "Some help string"),
                    }
  },

  "TestOptions_float" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of floats/doubles",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                   Type    Default     Req     Check   Help
                        "float_opt_def":    (  OptionType.Float,  5.0,       False,   None,  "Some help string"),
                        "float_req":        (  OptionType.Float,  None,      True,    None,  "Some help string"),
                        "float_opt":        (  OptionType.Float,  None,      False,   None,  "Some help string"),
                    }
  },

  "TestOptions_bool" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of bools",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                        "bool_opt_def":    (  OptionType.Bool,  False,     False,   None,  "Some help string"),
                        "bool_req":        (  OptionType.Bool,  None,      True,    None,  "Some help string"),
                        "bool_opt":        (  OptionType.Bool,  None,      False,   None,  "Some help string"),
                    }
  },

  "TestOptions_str" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of strings",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type    Default     Req     Check   Help
                        "str_opt_def":     (  OptionType.String,   "Hello",   False,   None,  "Some help string"),
                        "str_req":         (  OptionType.String,   None,      True,    None,  "Some help string"),
                        "str_opt":         (  OptionType.String,   None,      False,   None,  "Some help string"),
                    }
  },

  "TestOptions_listint" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of list of integers",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type         Default         Req     Check   Help
                        "listint_opt_def":  (  OptionType.ListInt,   [ 1, 2, 3 ],   False,   None,  "Some help string"),
                        "listint_req":      (  OptionType.ListInt,   None,          True,    None,  "Some help string"),
                        "listint_opt":      (  OptionType.ListInt,   None,          False,   None,  "Some help string"),
                    }
  },

  "TestOptions_listfloat" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of list of floats/doubles",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                   Type         Default                Req      Check   Help
                        "listfloat_opt_def":  (  OptionType.ListFloat,   [ 1.0, 2.0, 3.0 ],   False,   None,  "Some help string"),
                        "listfloat_req":      (  OptionType.ListFloat,   None,                True,    None,  "Some help string"),
                        "listfloat_opt":      (  OptionType.ListFloat,   None,                False,   None,  "Some help string"),
                    }
  },

  "TestOptions_listbool" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of list of bools",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type         Default                Req    Check   Help
                        "listbool_opt_def":  ( OptionType.ListBool, [ True, False, True ],  False, None, "Some help string"),
                        "listbool_req":      ( OptionType.ListBool,   None,                 True,  None, "Some help string"),
                        "listbool_opt":      ( OptionType.ListBool,   None,                 False, None, "Some help string"),
                    }
  },

  "TestOptions_liststr" :
  {
    "type"        : "python_module",
    "version"     : "0.1a",
    "description" : "Tests option parsing of list of strings",
    "authors"     : [],
    "refs"        : [],
    "options"     : {
                        # Key                 Type         Default                   Req     Check   Help
                        "liststr_opt_def":  ( OptionType.ListString,    [ "A", "Test", "List" ],  False, None, "Some help string"),
                        "liststr_req":      ( OptionType.ListString,    None,                     True,  None, "Some help string"),
                        "liststr_opt":      ( OptionType.ListString,    None,                     False, None, "Some help string"),
                    }
  },


}


