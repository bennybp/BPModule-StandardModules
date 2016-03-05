import bpmodule as bp

class TestOptions(bp.modulebase.Test_Base):
  def __init__(self, myid):
    super(TestOptions, self).__init__(self, myid)


  def RunTest_(self):
    bp.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def CallRunTest_(self, s):
    bp.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def Throw_(self):
    bp.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def CallThrow_(self, s):
    bp.output.Output("+++ In TestOptions. This is not meant to be run\n")

