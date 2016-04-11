import bpmodule as bp
from subprocess import call


def CallPsi4(EMethod,MethodName,order):
      Mol=EMethod.InitialWfn().system
      f=open(MethodName+".in","w")
      f.write("molecule{\n")
      f.write("units=bohr\n")
      f.write("no_reorient\n")
      f.write("no_com\n")
      NAtoms=0
      for atom in Mol:
        Symbol=atom.GetSymbol()
        if Symbol == "Gho" :
           if NAtoms%3==0:
              f.write("@O")
           else:
              f.write("@H")
            
        else:
           f.write(Symbol)
        f.write(" "+str(atom[0])+" "+str(atom[1])+" "+str(atom[2])+"\n")
        NAtoms+=1
      f.write("}\n")
      f.write("memory 10 gb\n")
      f.write("set{\n")
      f.write("freeze_core true\n")
      print(EMethod.Options().Get("BASIS_SET"))
      f.write("basis "+EMethod.Options().Get("BASIS_SET")+"\n")
      f.write("}\n")
      if order==0:
          f.write("energy(\""+MethodName+"\")\n")
      elif order==1:
          f.write("gradient(\""+MethodName+"\")\n")
      f.close()
      call(["/theoryfs2/ds/richard/SrcFiles/psi4/objdir/bin/psi4",MethodName+".in"])
      f=open(MethodName+".out","r")
      egy=[]
      for line in f:
         if order==0:
            if MethodName=="SCF":
               if line.split()[0:3]==["@DF-RHF","Final","Energy:"]:
                  egy.append(float(line.split()[3]))
                  break
            if MethodName=="MP2":
               if line.split()[0:3]==["Total","Energy","="]:
                  egy.append(float(line.split()[3]))
                  break
         elif order==1 and line.split()[0:2]==["-Total","Gradient:"]:
             line=next(f)
             line=next(f)
             line=next(f)
             for x in range(0,NAtoms):
                 for y in range(1,4):
                    egy.append(float(line.split()[y]))
                 line=next(f)
             break
      f.close()
      return egy
     
