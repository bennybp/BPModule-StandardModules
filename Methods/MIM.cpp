#include <stack>
#include "Methods/MIM/MIM.hpp" //Note that MIM.hpp is not in a folder yet b/c I can't figure out how to make folders on GitHub

using bpmodule::modulemanager::ModuleManager;
using bpmodule::system::System;
using bpmodule::system::SystemMap;

typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;
typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_;
typedef std::vector<double> Return_t;
typedef std::map<std::string,Return_t> DerivMap;

namespace bpmethod{

//TODO: pull from MBE.cpp
class Task{
   private:
      ModuleManager& MM_;
      std::string Key_;
      unsigned long ID_;
      const System& Sys_;
   public:
      Task(ModuleManager& MM, const std::string& Key, unsigned long ID, const System& Sys):
         MM_(MM),Key_(Key),ID_(ID),Sys_(Sys){}
      Return_t operator()(size_t Order)const{
         EMethod_t DaMethod=MM_.GetModule<EnergyMethod>(Key_,ID_);
         DaMethod->Wfn().system.Set(Sys_);
         return DaMethod->Deriv(Order);
      }
};

/*
 *   Mapping atoms from subsystems to the supersystem, always fun.  Here's the plan.  For the Order-th
 *   derivative we have to expand a rank "Order" tensor into another rank "Order" tensor.  We assume that
 *   the order of elements is x, y, z for atom1; x,y,z for atom2; etc.  Where atom1,atom2, etc. are relative
 *   to the order of atoms specified in the incoming subsystem.  Technically the derivatives have symmetry,
 *   but we ignore that for now.  We do the actual filling by recursion.  At each level we loop over atoms
 *   and components and pass the result to the next level.  Once we have gone down Order levels, we have a fully
 *   specified index and simply evaluate the index.  We can avoid a little-bit of recursion (and possibly pick up
 *   a tiny bit of vectorization if we unroll the last iteration.  We thus will have two offsets.  The first
 *   offset is the offset of the Order-1 indices.  The second offset is the order-th index.  Note that each order
 *   is sorta two orders because we get an offset from the components too.
 *   Basically our index looks like:
 *   Idx=\left[\sum_{i=0}^{Order-1} Atom[i]*(3*NAtoms)^{Order-1-i}+Comp[i]\right]+AtomJ*3+CompJ;
 */
void FillDeriv(Result_t& Result, 
               const Result_t& SubResult,
               double Coeff,
               const System& Sys, 
               std::stack<Atom> Idx,
               std::stack<size_t> Comp,
               const std::map<Atom,size_t>& SuperAtomMap,
               const std::map<Atom,size_t>& SubAtomMap,
               size_t Order){
   //Handle energy
   if(Order==0)Result[0]+=Coeff*SubResult[0];
   else if(Idx.size()==(Order-1)){//Gradient lands here, Hessian and higer, land here on recursion
      size_t SuperOff=0,SubOff=0;
      for(size_t i=0;i<Order-1;++i){//First offset
         size_t SuperStride=1,SubStride=1;
         for(size_t j=i;j<Order-1;++j){
            SuperStride*=3*SuperAtomMap.size();
            SubStride*=3*SubAtomMap.size();
         }
         SuperOff+=(SuperAtomMap[Idx[i]]+Comp[i])*SuperStride;
         SubOff+=(SubAtomMap[Idx[i]]+Comp[i])*SubStride;
      }
      for(const Atom& AtomI: Sys){//Unrolled loop
         //Second offsets
         size_t SuperOff2=SuperAtomMap[AtomI]*3;
         size_t SubOff2=SubAtomMap[AtomI]*3;
         for(size_t i=0;i<3;++i)//Actual filling
            Result[SuperOff+SuperOff2+i]+=Coeff*SubResult[SubOff+SubOff2+i];
      }
   }
   else{//Hessian and higher land here
      for(const Atom& AtomI : Sys){
         Idx.push(AtomI);
         for(size_t i=0;i<3;++i){
            Comp.push(i);
            FillDeriv(Result,SubResult,Coeff,Sys,Idx,Comp,SuperAtomMap,SubAtomMap,Order);
            Comp.pop();
         }
         Idx.pop();
      }
   }
   
}

Return_t MIM::DerivImpl(size_t Order)const{
   //Get the system and compute the number of degrees of freedom for the result
   const System& Mol=*Wfn.system;
   size_t DoF=1;
   for(size_t i=0;i<Order;++i)DoF*=3*Mol.NAtoms();
   
   //Establish an atom order
   std::map<Atom,size_t> AtomMap;
   size_t counter=0
   for(const Atom& AtomI : Mol)AtomMap[AtomI]=counter++;
   
   //Get the options
   const OptionMap& DaOptions=Options();
   std::vector<std::string> MethodNames=DaOptions.Get<std::vector<std::string>>("METHODS");
   Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>("FRAG");
   SystemMap Systems=Fragger->Fragmentize(Mol);
   std::vector<double> Coeffs=DaOptions.Get<std::vector<double>>("WEIGHTS");
   
   //For the time-being the user is required to give us a coefficient for each task
   size_t NTasks=Coeffs.size();
   
   //True if we are using the same method for all systems
   bool SameMethod=MethodNames.size()==1;
   //True if we are using the same system for all methods
   bool SameSystem=Systems.size()==1;
   
   const Communicator& ParentComm=bpmodule::parallel::GetEnv().Comm();
   //For the time-being we assume that we are only using as many processes as methods
   Communicator NewComm=ParentComm.Split(ParentComm.NThreads(),1,Systems.size());
   TaskResults<Return_t> Results(NewComm);
   
   //Use two loops so all tasks get queued before we start asking for results
   for(size_t TaskI=0; TaskI<NTasks; ++TaskI){
      Results.push_back(
         NewComm.AddTask(
            Task(MManager(),
                 MethodNames[SameMethods?0:TaskI],
                 ID(),
                 Systems[SameSystem?0:TaskI].second;
            ),
         Order);
   }
   
   //Will be the derivatives per system
   DerivMap Derivs;
   //Final, total derivative
   Result_t TotalDeriv(DoF,0.0);
   //Buffers for the recursion
   std::stack<Atom> Idx;
   std::stack<size_t> Comp;
   
   for(size_t TaskI=0;TaskI<NTasks;++TaskI){
      const System& SubSys=Systems[SameSystem?0:TaskI].first;
      Derivs[SubSys]=Results[TaskI];
      std::map<Atom,size_t> SubAtomMap;
      counter=0;
      for(const Atom& AtomI: SubSys)SubAtomMap[SubSys]=counter++;
      FillDeriv(TotalDeriv,Results[TaskI],Coffs[TaskI],SubSys,Idx,Comp,AtomMap,SubAtomMap,Order);
   }
   return TotalDeriv;
}

}//End namespace

