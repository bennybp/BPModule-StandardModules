import pulsar as psr

class TestOptions(psr.Test_Base):
  def __init__(self, myid):
    super(TestOptions, self).__init__(self, myid)


  def run_test_(self):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def call_run_test_(self, s):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def Throw_(self):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

  def call_throw_(self, s):
    psr.output.Output("+++ In TestOptions. This is not meant to be run\n")

