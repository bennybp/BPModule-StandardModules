import pulsar as psr

class TestPyModule1(psr.modulebase.Test_Base):
  def __init__(self, myid):
    super(TestPyModule1, self).__init__(myid)


  def run_test_(self):
    self.out.output("+++ In TestPyModule1: run_test. Info: ({}) {} {} v{}\n".format(self.id(), self.key(), self.name(), self.version()))

    self.out.output("   Cache entries: {}\n".format(self.cache().size()))

    self.cache().print(self.out)

    self.out.output("   double_opt_def:    {}\n".format(self.options().get("double_opt_def")))
    self.out.output("      int_opt_def:    {}\n".format(self.options().get("int_opt_def")))
    self.out.output("     bool_opt_def:    {}\n".format(self.options().get("bool_opt_def")))
    self.out.output("      str_opt_def:    {}\n".format(self.options().get("str_opt_def")))
    self.out.output("\n")
    if self.options().has("double_opt"):
        self.out.output("       double_opt:    {}\n".format(self.options().get("double_opt")))
    if self.options().has("int_opt"):
        self.out.output("          int_opt:    {}\n".format(self.options().get("int_opt")))
    if self.options().has("bool_opt"):
        self.out.output("         bool_opt:    {}\n".format(self.options().get("bool_opt")))
    if self.options().has("str_opt"):
        self.out.output("          str_opt:    {}\n".format(self.options().get("str_opt")))


    # cache something
    cp = int(psr.datastore.CacheData.CheckpointLocal)

    self.cache().set( "Element 1", "Something in the python cache", cp)
    self.cache().set( "Element 2", 42, cp)
    self.cache().set( "Element 3", 42.0, cp)
    self.cache().set( "Element 4", [ 1, 2, 3, 4 ], cp)

  def call_run_test_(self, other):
    self.out.output("+++ In TestPyModule1: call_run_test with {}\n".format(other))

    tb = self.create_child(other)
    self.out.output("  + Obtained scoped module ID {}\n".format(tb.id()))
    tb.run_test()
    self.out.output("  + Finished with scoped module {}. Deleting automatically\n".format(tb.id()))

    self.out.output("+++Done\n")


  def call_run_test2_(self, other1, other2):
    self.out.output("+++ In TestPyModule1: call_run_test with {} {}\n".format(other1, other2))

    tb = self.create_child(other1)
    self.out.output("  + Obtained scoped module ID {}\n".format(tb.id()))
    tb.call_run_test(other2)
    self.out.output("  + Finished with scoped module {}. Deleting automatically\n".format(tb.id()))

    self.out.output("+++Done\n")

  def test_throw_(self):
    psr.output.GlobalWarning("+++ In TestPyModule1: Throwing an exception!\n")
    raise psr.exception.GeneralException("Here in py", "Key", "Some Data")
    #self.Throw("This is a test exception from python")


  def call_throw_(self, other):
    self.out.output("+++ In TestPyModule1: call_run_test with {}\n".format(other))

    tb = self.create_child(other)
    self.out.output("  + Obtained scoped module ID {}\n".format(tb.id()))
    tb.test_throw()

    # shouldn't be run?
    self.out.output("+++Done\n")
