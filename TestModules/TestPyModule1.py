import bpmodule as bp

class TestPyModule1(bp.modulebase.Test_Base):
  def __init__(self, myid):
    super(TestPyModule1, self).__init__(myid)


  def RunTest_(self):
    bp.output.GlobalOutput("+++ In TestPyModule1: RunTest. Info: ({}) {} {} v{}\n".format(self.ID(), self.Key(), self.Name(), self.Version()))

    bp.output.GlobalOutput("    Wavefunction: {}\n".format(self.InitialWfn().UniqueString()))
    bp.output.GlobalOutput("   Cache entries: {}\n".format(self.Cache().Size()))

    for it in self.Cache().GetKeys():
        bp.output.GlobalOutput("                  > {}\n".format(it))

    bp.output.GlobalOutput("   double_opt_def:    {}\n".format(self.Options().Get("double_opt_def")))
    bp.output.GlobalOutput("      int_opt_def:    {}\n".format(self.Options().Get("int_opt_def")))
    bp.output.GlobalOutput("     bool_opt_def:    {}\n".format(self.Options().Get("bool_opt_def")))
    bp.output.GlobalOutput("      str_opt_def:    {}\n".format(self.Options().Get("str_opt_def")))
    bp.output.GlobalOutput("\n")
    if self.Options().Has("double_opt"):
        bp.output.GlobalOutput("       double_opt:    {}\n".format(self.Options().Get("double_opt")))
    if self.Options().Has("int_opt"):
        bp.output.GlobalOutput("          int_opt:    {}\n".format(self.Options().Get("int_opt")))
    if self.Options().Has("bool_opt"):
        bp.output.GlobalOutput("         bool_opt:    {}\n".format(self.Options().Get("bool_opt")))
    if self.Options().Has("str_opt"):
        bp.output.GlobalOutput("          str_opt:    {}\n".format(self.Options().Get("str_opt")))


    # cache something
    self.Cache().Set( "Element 1", "Something in the python cache")
    self.Cache().Set( "Element 2", 42)
    self.Cache().Set( "Element 3", 42.0)
    self.Cache().Set( "Element 4", [ 1, 2, 3, 4 ])

  def CallRunTest_(self, other):
    bp.output.GlobalOutput("+++ In TestPyModule1: CallRunTest with {}\n".format(other))

    tb = self.CreateChildModule(other)
    bp.output.GlobalOutput("  + Obtained scoped module ID {}\n".format(tb.ID()))
    tb.RunTest()
    bp.output.GlobalOutput("  + Finished with scoped module {}. Deleting automatically\n".format(tb.ID()))

    bp.output.GlobalOutput("+++Done\n")


  def CallRunTest2_(self, other1, other2):
    bp.output.GlobalOutput("+++ In TestPyModule1: CallRunTest with {} {}\n".format(other1, other2))

    tb = self.CreateChildModule(other1)
    bp.output.GlobalOutput("  + Obtained scoped module ID {}\n".format(tb.ID()))
    tb.CallRunTest(other2)
    bp.output.GlobalOutput("  + Finished with scoped module {}. Deleting automatically\n".format(tb.ID()))

    bp.output.GlobalOutput("+++Done\n")

  def TestThrow_(self):
    bp.output.GlobalWarning("+++ In TestPyModule1: Throwing an exception!\n")
    raise bp.exception.GeneralException("Here in py", "Key", "Some Data")
    #self.Throw("This is a test exception from python")


  def CallThrow_(self, other):
    bp.output.GlobalOutput("+++ In TestPyModule1: CallRunTest with {}\n".format(other))

    tb = self.CreateChildModule(other)
    bp.output.GlobalOutput("  + Obtained scoped module ID {}\n".format(tb.ID()))
    tb.TestThrow()

    # shouldn't be run?
    bp.output.GlobalOutput("+++Done\n")
