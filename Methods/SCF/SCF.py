import bpmodule as bp
from subprocess import call

class SCF(bp.modulebase.EnergyMethod):
  def __init__(self, myid):
    super(SCF, self).__init__(myid)
  
  def Deriv_(self,order):
      Mol=self.Wfn().system.Get()
      f=open("MBE.in","w")
      f.write("molecule{\n")
      f.write("units=bohr\n")
      for atom in Mol:
          f.write(str(atom.GetSymbol())+" "+str(atom[0])+" "+str(atom[1])+
                    " "+str(atom[2])+"\n")
      f.write("}\n")
      f.write("memory 10 gb\n")
      f.write("set{\n")
      f.write("basis sto-3g\n")
      f.write("}\n")
      f.write("energy(\"SCF\")\n")
      f.close()
      call(["/theoryfs2/ds/richard/SrcFiles/psi4/objdir/bin/psi4","MBE.in"])
      f=open("MBE.out","r")
      egy=0.0
      for line in f:
         if line.split()[0:3]==["@DF-RHF","Final","Energy:"]:
             egy=line.split()[3]
      f.close()
      return [float(egy)]
     
