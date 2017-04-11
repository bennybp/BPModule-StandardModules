import pulsar as psr
import numpy as np
from scipy.optimize import *

class geom_functor:
    def __init__(self,wfn,mod):
        self.wfn_=wfn
        self.mod_=mod
    def make_wfn(self,geom):
        new_uv=psr.AtomSetUniverse()
        for a,carts in zip(self.wfn_.system,zip(*[iter(geom)]*3)):
            temp_atom=psr.Atom(a)
            for i in range(3):
                temp_atom[i]=carts[i]
            new_uv.insert(temp_atom)
        self.wfn_.system=psr.System(new_uv,True)
        return self.wfn_
    def grad(self,geom):
        new_wfn,grad=self.mod_.deriv(1,self.make_wfn(geom))
        self.wfn_=new_wfn
        print("Max Absolute Gradient Component :",max(abs(x) for x in grad))
        return np.array(grad)
    def __call__(self, geom):
        temp,egy=self.mod_.deriv(0,self.make_wfn(geom))
        print("New Energy:",egy[0])
        return egy[0]

class GeometryOptimizer(psr.EnergyMethod):
    """This class is a thin wrapper around the SciPy optimizer
    """

    def __init__(self,myid):
        super(GeometryOptimizer,self).__init__(myid)

    def deriv_(self,order,wfn):
        mod=self.create_child_from_option("METHOD_KEY")
        fxn=geom_functor(wfn,mod)
        x0=[]
        for a in wfn.system:
            for i in range(3):
                x0.append(a[i])
        GTol=self.options().get("MAX_GRAD")
        meth=self.options().get("TYPE")
        opt={"disp":True}
        result=minimize(fxn,x0,method=meth,jac=fxn.grad,tol=GTol,options=opt)
        return mod.deriv(order,fxn.make_wfn(result.x))
