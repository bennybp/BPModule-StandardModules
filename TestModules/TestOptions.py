import pulsar as psr

class TestOptions(psr.modulebase.Test_Base):
  def __init__(self, myid):
    super(TestOptions, self).__init__(self, myid)


  def RunTest_(self):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def CallRunTest_(self, s):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def Throw_(self):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def CallThrow_(self, s):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

