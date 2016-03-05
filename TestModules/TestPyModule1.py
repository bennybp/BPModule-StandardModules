import bpmodule as bp

class TestPyModule1(bp.modulebase.Test_Base):
  def __init__(self, myid):
    super(TestPyModule1, self).__init__(myid)


  def RunTest_(self):
    bp.output.Output("+++ In TestPyModule1: RunTest. Info: (%1%) %2% %3% v%4%\n", self.ID(), self.Key(), self.Name(), self.Version());

    bp.output.Output("    Wavefunction: %1% \n", self.Wfn().UniqueString())
    bp.output.Output("   Cache entries: %1%\n", self.Cache().Size());

    for it in self.Cache().GetKeys():
        bp.output.Output("                  > %1%\n", it);

    bp.output.Output("   double_opt_def:    %1%\n", self.Options().Get("double_opt_def"));
    bp.output.Output("      int_opt_def:    %1%\n", self.Options().Get("int_opt_def"));
    bp.output.Output("     bool_opt_def:    %1%\n", self.Options().Get("bool_opt_def"));
    bp.output.Output("      str_opt_def:    %1%\n", self.Options().Get("str_opt_def"));
    bp.output.Output("\n");
    if self.Options().Has("double_opt"):
        bp.output.Output("       double_opt:    %1%\n", self.Options().Get("double_opt"));
    if self.Options().Has("int_opt"):
        bp.output.Output("          int_opt:    %1%\n", self.Options().Get("int_opt"));
    if self.Options().Has("bool_opt"):
        bp.output.Output("         bool_opt:    %1%\n", self.Options().Get("bool_opt"));
    if self.Options().Has("str_opt"):
        bp.output.Output("          str_opt:    %1%\n", self.Options().Get("str_opt"));


    # cache something
    self.Cache().Set( "Element 1", "Something in the python cache")
    self.Cache().Set( "Element 2", 42)
    self.Cache().Set( "Element 3", 42.0)
    self.Cache().Set( "Element 4", [ 1, 2, 3, 4 ])

  def CallRunTest_(self, s):
    bp.output.Output("+++ In TestPyModule1: CallRunTest with %1%\n", s)

    tb = self.CreateChildModule(s)
    bp.output.Output("  + Obtained scoped module ID %1%\n", tb.ID())
    tb.RunTest()
    bp.output.Output("  + Finished with scoped module %1%. Deleting automatically\n", tb.ID())

    bp.output.Output("+++Done\n");

  def TestThrow_(self):
    bp.output.Warning("+++ In TestPyModule1: Throwing an exception!\n");
    raise bp.exception.GeneralException("Here in py", "Key", "Some Data")
    #self.Throw("This is a test exception from python")


  def CallThrow_(self, s):
    bp.output.Output("+++ In TestPyModule1: CallRunTest with %1%\n", s)

    tb = self.CreateChildModule(s)
    bp.output.Output("  + Obtained scoped module ID %1%\n", tb.ID())
    tb.TestThrow()

    # shouldn't be run?
    bp.output.Output("+++Done\n");
