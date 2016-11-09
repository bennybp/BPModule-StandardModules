import pulsar as psr

def Grad(NDeriv,DDeriv,Egy):
   G=[]
   for i in range(0,len(NDeriv[1])):
      G.append((2*NDeriv[0][0]*NDeriv[1][i]-Egy*DDeriv[1][i])/DDeriv[0][0])
   return G

def Hess(NDeriv,DDeriv,Egy,G):
   H=[]
   Stride=len(NDeriv[1])
   for i in range(0,Stride):
      for j in range(0,Stride):
          H.append((2*NDeriv[1][i]*NDeriv[1][j]+
                   2*NDeriv[0][0]*NDeriv[2][i*stride+j]-
                   Egy*DDeriv[2][i*stride+j]-
                   2*DDeriv[1][i]*G[j])/DDeriv[0][0])
   return H



class FellerCBS(psr.EnergyMethod):
  """This module implements the three-point Feller extrapolation

  In the Feller extrapolation scheme, the energy for a basis set of cardinal 
  number \f$X\f$ is given by:
  \f[
    E_X=E_{CBS}+\alpha\exp(-\beta X).
  \f]
  For three basis sets of cardinal numbers: \f$X\f$,\f$X+1\f$, and \f$X+2\f$,
  the equations can be rearranged to give:
  \f[
     E_{CBS}=E_X-\frac{\left(E_X-E_{X+1}\right)^2}{E_X-2E_{X+1}+E_{X+2}}=
E_X-\frac{E^2_N}{E_D}.
  \f]
  where we have defined the following:
  \f{align*}{
    E_{N}\equiv&E_X-E_{X+1}\\
    E_{D}\equiv& E_X-2E_{X+1}+E_{X+2}\\
    \epsilon\equiv& E_X-E_{CBS}.
  \f}
  For derivatives, it is easier if we rearrange this so the linear terms are on
  one side:
  \f[
     \epsilon=E_X-E_{CBS}=\frac{E_N^2}{E_D}.
  \f]

  The gradient of this is given by:
    \f{align*}{
     \bigtriangledown \epsilon=&
     \frac{2E_N \bigtriangledown E_N}{E_D}-
      \frac{E^2_N\bigtriangledown E_D}{E^2_D}\\
      =&     \frac{2E_N \bigtriangledown E_N}{E_D}-
      \frac{\epsilon\bigtriangledown E_D}{E_D}\\
      =&\frac{2E_N \bigtriangledown E_N-\epsilon\bigtriangledown E_D}{E_D}\\
  \f}
  The Hessian is (multiplication between gradients is tensor product):
  \f{align*}{
     \epsilon^{(2)}=&\frac{2 E^2_N \left(E^{(1)}_D\right)^2}{E_D^3}-
                  \frac{4E_N E^{(1)}_N E^{(1)}_D}{E_D^2}+
                  \frac{2 \left(E_N^{(1)}\right)^2}{E_D}-
                  \frac{E_N^2 E^{(2)}_D}{E_D^2}+
                  \frac{2 E_N E_N^{(2)}}{E_D}\\
                  =&\frac{2 \epsilon \left(E^{(1)}_D\right)^2}{E_D^2}-
                  \frac{4E_N E^{(1)}_N E^{(1)}_D}{E_D^2}+
                  \frac{2 \left(E_N^{(1)}\right)^2}{E_D}-
                  \frac{\epsilon E^{(2)}_D}{E_D}+
                  \frac{2 E_N E_N^{(2)}}{E_D}\\
                  =&\frac{2\left(\epsilon E^{(1)}_D-
                          2 E_N E^{(1)}_N\right)E^{(1)}_D}{E_D^2}+
                  \frac{2 \left(E_N^{(1)}\right)^2+2 E_N E_N^{(2)}-
                        \epsilon E^{(2)}_D}{E_D}\\
                  =&\frac{2 \left(E_N^{(1)}\right)^2+2 E_N E_N^{(2)}-
                          2\epsilon^{(1)}E^{(1)}_D-
                        \epsilon E^{(2)}_D}{E_D}
  \f}

  """
  def __init__(self, myid):
    super(FellerCBS, self).__init__(myid)

  def deriv_(self,order,wfn):
      #Get derivatives of numerator, denominator, and smallest basis set
      Derivs[{},{},{}]
      Cs=[[1.0],[1.0,-1.0],[1.0,-2.0,1.0]]
      #Do Denom first as num calcs are subset
      for i in range(2,-1):
         MIM=self.GetSubModule(self.options().get("MIM_KEY"))
         MIM.options().change("WEIGHTS",Cs[i])
         MIM.options().change("METHODS",[Options.get("METHOD")])       
         Deriv[i][order]=MIM.deriv(order)
         for j in range(0,order-1):
            Deriv[j][order]=MIM.deriv(j)
      Epsilon=Deriv[1][0][0]**2/Deriv[2][0][0]
      if order==0 :
         return [Deriv[0][0][0]-Epsilon]
      G=Grad(Deriv[1],Deriv[2],Epsilon)
      if order==1:
         return [Deriv[0][1][i]-G[i] for i in G]
      H=Hess(Deriv[1],Deriv[2],Epsilon,G)
      return (Wfn,[Deriv[0][2][i]-H[i] for i in H])
         
      

      
