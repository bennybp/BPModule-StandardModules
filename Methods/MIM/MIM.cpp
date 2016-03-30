#include <vector>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <bpmodule/system/System.hpp>
#include <bpmodule/exception/Exceptions.hpp>
#include <bpmodule/modulebase/SystemFragmenter.hpp>
#include <bpmodule/datastore/OptionMap.hpp>
#include <bpmodule/parallel/InitFinalize.hpp>
#include <LibTaskForce.hpp>
#include <bpmodule/output/Table.hpp>
#include "Methods/MIM/MIM.hpp" 

using bpmodule::modulemanager::ModuleManager;
using bpmodule::system::System;
using bpmodule::system::SystemMap;
using bpmodule::system::Atom;
using bpmodule::modulemanager::ModuleManager;
using bpmodule::modulebase::SystemFragmenter;
using bpmodule::modulebase::EnergyMethod;
using bpmodule::datastore::OptionMap;
using LibTaskForce::Communicator;
using LibTaskForce::TaskResults;
using std::vector;
using std::map;
using std::unordered_map;
using std::string;


typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;
typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;
typedef vector<double> Return_t;
typedef map<string,Return_t> DerivMap;
typedef unsigned long ULI;
typedef unordered_map<Atom,size_t> AtomMap_t;

namespace bpmethods{

class Task{
   private:
      ModuleManager& MM_;
      string Key_;
      ULI ID_;
      const System& Sys_;
   public:
      Task(ModuleManager& MM, const string& Key, ULI ID, const System& Sys):
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
void FillDeriv(Return_t& Result, 
               const Return_t& SubResult,
               double Coeff,
               const System& Sys, 
               const AtomMap_t& SuperAtomMap,
               const AtomMap_t& SubAtomMap,
               size_t Order,
               vector<Atom> Idx=vector<Atom>(),
               vector<size_t> Comp=vector<size_t>()){
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
         SuperOff+=(SuperAtomMap.at(Idx[i])+Comp[i])*SuperStride;
         SubOff+=(SubAtomMap.at(Idx[i])+Comp[i])*SubStride;
      }
      for(const Atom& AtomI: Sys){//Unrolled loop
         //Second offsets
         size_t SuperOff2=SuperAtomMap.at(AtomI)*3;
         size_t SubOff2=SubAtomMap.at(AtomI)*3;
         for(size_t i=0;i<3;++i)//Actual filling
            Result[SuperOff+SuperOff2+i]+=Coeff*SubResult[SubOff+SubOff2+i];
      }
   }
   else{//Hessian and higher land here
      for(const Atom& AtomI : Sys){
         Idx.push_back(AtomI);
         for(size_t i=0;i<3;++i){
            Comp.push_back(i);
            FillDeriv(Result,SubResult,Coeff,Sys,
                      SuperAtomMap,SubAtomMap,Order,Idx,Comp);
            Comp.pop_back();
         }
         Idx.pop_back();
      }
   }
   
}

void PrintTable(const vector<string>& Rows,const DerivMap& Derivs);

AtomMap_t MapAtoms(const System& Mol){
   unordered_map<Atom,size_t> AtomMap;
   size_t counter=0;
   for(const Atom& AtomI : Mol)AtomMap[AtomI]=counter++;
   return AtomMap;
}

Return_t MIM::DerivImpl(size_t Order)const{
   //Get the system and compute the number of degrees of freedom for the result
   const System& Mol=*Wfn().system;
   size_t DoF=1;
   for(size_t i=0;i<Order;++i)DoF*=3*Mol.NAtoms();
   
   //Establish an atom order
   AtomMap_t AtomMap=MapAtoms(Mol);
   
   //Get the options
   const OptionMap& DaOptions=Options();
   vector<string> MethodNames=DaOptions.Get<vector<string>>("METHODS");
   string Fragmentizer=DaOptions.Get<string>("FRAGMENTIZER");
   Return_t Coeffs=DaOptions.Get<Return_t>("WEIGHTS");
   
   //Get the systems
   Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>(Fragmentizer);
   SystemMap Systems=Fragger->Fragmentize(Mol);
   
   //For the time-being the user is required to give us a coefficient for each task
   size_t NTasks=Coeffs.size();
   
   //True if we are using the same method for all systems
   bool SameMethod=MethodNames.size()==1;
   //True if we are using the same system for all methods
   bool SameSystem=Systems.size()==1;
   
   if(SameSystem && SameMethod && (NTasks>1))
       throw bpmodule::exception::GeneralException(
               "Minimally, either the number of systems or the number of "
               "methods must equal the number of coefficients");
   
   //Set-up parallel and our buffer
   const Communicator& ParentComm=bpmodule::parallel::GetEnv().Comm();
   Communicator NewComm=ParentComm.Split(ParentComm.NThreads(),1,
                                         std::min(ParentComm.NProcs(),NTasks));
   TaskResults<Return_t> Results(NewComm);
   
   //Use two loops. Gets all tasks queued before we start asking for results
   SystemMap::const_iterator SysI=Systems.begin();
   for(size_t TaskI=0; TaskI<NTasks; ++TaskI){
      Results.push_back(
         NewComm.AddTask(
            Task(MManager(),MethodNames[SameMethod?0:TaskI],ID(),SysI->second),
            Order
         )
      );
      if(!SameSystem)++SysI;
   }
   
   
   //Loop two, gettin results
   DerivMap Derivs; //Will be the derivatives per system
   Return_t TotalDeriv(DoF,0.0);//Final, total derivative
   
   //Iterator to our first system
   SysI=Systems.begin();
   vector<string> RowTitles(NTasks);
   for(size_t TaskI=0;TaskI<NTasks;++TaskI){
      const System& SubSys=SysI->second;
      Derivs[SysI->first]=Results[TaskI];
      RowTitles[TaskI]=SysI->first;
      FillDeriv(TotalDeriv,Results[TaskI],Coeffs[TaskI],SubSys,
                AtomMap,MapAtoms(SubSys),Order);
      if(!SameSystem)++SysI;
   }
   if(Order==0)PrintTable(RowTitles,Derivs);
   return TotalDeriv;
}

void PrintTable(const vector<string>& Rows,const DerivMap& Derivs){
   const size_t NCols=2,NRows=Rows.size()+1;
   bpmodule::output::Table ResultTable(NRows,NCols);
   std::array<string,NCols> ColTitles={"System","Energy"};
   ResultTable.SetHBorder(0,'*');
   ResultTable.SetHBorder(1,'-');
   ResultTable.SetHBorder(NRows,'*');
   ResultTable.SetVBorder(1,'|');
   ResultTable.FillRow(ColTitles,0,0,NCols);
   ResultTable.FillCol(Rows,0,1,NRows);
   for(size_t i=0;i<Rows.size();++i)
         ResultTable.FillRow(&Derivs.at(Rows[i])[0],i+1,1,NCols);
   //bpmodule::output::GetGlobalOut()<<ResultTable<<std::endl;
}


}//End namespace

