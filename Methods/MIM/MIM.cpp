#include <unordered_map>
#include <limits>
#include <algorithm>
#include <memory>

#include <pulsar/exception/Exceptions.hpp>
#include <pulsar/parallel/InitFinalize.hpp>
#include <LibTaskForce.hpp>
#include <pulsar/output/Table.hpp>
#include <pulsar/output/GlobalOutput.hpp>
#include "Methods/MIM/MIM.hpp" 
#include "Methods/MBE/MBECommon.hpp"
#include "Methods/MBE/MBEUtils.hpp"

using pulsar::modulemanager::ModuleManager;
using LibTaskForce::Communicator;
using LibTaskForce::TaskResults;
using std::vector;
using std::map;
using std::unordered_map;
using std::string;

typedef vector<double> Return_t;
typedef map<string,Return_t> DerivMap;
typedef unsigned long ULI;
typedef unordered_map<Atom,size_t> AtomMap_t;

namespace pulsarmethods{

class Task{
   private:

      ModuleManager& MM_;
      string Key_;
      ULI ID_;
      const System& Sys_;
      size_t TaskNum_;
   public:
      Task(ModuleManager& MM, const string& Key, ULI ID, const System& Sys,
           size_t TaskNum=0):
         MM_(MM),Key_(Key),ID_(ID),Sys_(Sys),TaskNum_(TaskNum){}
      Return_t operator()(size_t Order)const{ 
        EMethod_t DaMethod=MM_.GetModule<EnergyMethod>(Key_,ID_);
        DaMethod->InitialWfn().GetSystem()=std::make_shared<System>(Sys_);
        return DaMethod->Deriv(Order);
      }
};



void PrintEgyTable(const vector<string>& Rows,const DerivMap& Derivs);
void PrintGradTable(const vector<string>& Rows,const DerivMap& Derivs,
                    const SystemMap& Systems);

AtomMap_t MapAtoms(const System& Mol){
   AtomMap_t AtomMap;
   size_t counter=0;
   for(const Atom& AtomI : Mol)AtomMap[AtomI]=counter++;
   return AtomMap;
}

Return_t MIM::DerivImpl(size_t Order)const{
   //Get the system and compute the number of degrees of freedom for the result
   const System& Mol=*InitialWfn().GetSystem();
   size_t DoF=1;
   for(size_t i=0;i<Order;++i)DoF*=3*Mol.Size();
    
   //Establish an atom order
   AtomMap_t AtomMap=MapAtoms(Mol);
   
   const OptionMap& DaOptions=Options();
   vector<string> MethodNames=DaOptions.Get<vector<string>>("METHODS");
   Return_t Coeffs=DaOptions.Get<Return_t>("WEIGHTS");
   

   //Get the subsystems
   Fragmenter_t Fragger=CreateChild<SystemFragmenter>(
           DaOptions.Get<string>("FRAGMENTIZER"));
   SystemMap Systems=Fragger->Fragmentize(Mol);
     
   //For the time-being the user is required to give us a coefficient for 
   //each task
   size_t NTasks=Coeffs.size();
   
   //True if we are using the same method for all systems
   bool SameMethod=MethodNames.size()==1;
   //True if we are using the same system for all methods
   bool SameSystem=Systems.size()==1;
   
   //TODO: move check to options
   if(SameSystem && SameMethod && (NTasks>1))
       throw pulsar::exception::GeneralException(
               "Minimally, either the number of systems "
               " or the number of methods must equal the number of coefficients");
   
   //Set-up parallel and our buffer
   const Communicator& ParentComm=pulsar::parallel::GetEnv().Comm();
   Communicator NewComm=ParentComm.Split(ParentComm.NThreads(),1,
                                         std::min(ParentComm.NProcs(),NTasks));
   TaskResults<Return_t> Results(NewComm);
   
   //Use two loops. Gets all tasks queued before we start asking for results
   SystemMap::const_iterator SysI=Systems.begin();
   for(size_t TaskI=0; TaskI<NTasks; ++TaskI){
      Results.push_back(
         NewComm.AddTask(
            Task(MManager(),
                MethodNames[SameMethod?0:TaskI],
                ID(),
                SysI->second
                //,NewComm.NTasks()
            ),
            Order
         )
      );
      if(!SameSystem)++SysI;
   }
   
   
   //Loop two, gettin results
   DerivMap Derivs; //Will be the derivatives per system
   Return_t TotalDeriv(DoF,0.0);//Final, total derivative
   
   SysI=Systems.begin();
   vector<string> RowTitles(NTasks);
   for(size_t TaskI=0;TaskI<NTasks;++TaskI){
      const System& SubSys=SysI->second;
      std::stringstream ss;
      ss<<SysI->first;
      ss<<" ["<<MethodNames[SameMethod?0:TaskI]<<"]";
      Derivs[ss.str()]=Results[TaskI];
      RowTitles[TaskI]=ss.str();
      FillDeriv(TotalDeriv,Results[TaskI],Coeffs[TaskI],SubSys,
                AtomMap,MapAtoms(SubSys),Order);
      if(!SameSystem)++SysI;
   }
   if(Order==0)PrintEgyTable(RowTitles,Derivs);
   else if(Order==1)PrintGradTable(RowTitles,Derivs,Systems);
   return TotalDeriv;
}

void PrintEgyTable(const vector<string>& Rows,const DerivMap& Derivs){
   const size_t NCols=2,NRows=Rows.size()+1;
   pulsar::output::Table ResultTable(NRows,NCols);
   std::array<string,NCols> ColTitles={"System [Method Key]","Energy (a.u.)"};
   ResultTable.SetHBorder(0,'*');
   ResultTable.SetHBorder(1,'-');
   ResultTable.SetHBorder(NRows,'*');
   ResultTable.SetVBorder(1,'|');
   ResultTable.FillRow(ColTitles,0,0,NCols);
   ResultTable.FillCol(Rows,0,1,NRows);
   for(size_t i=0;i<Rows.size();++i)
         ResultTable.FillRow(&Derivs.at(Rows[i])[0],i+1,1,NCols);
   pulsar::output::GetGlobalOut()<<ResultTable<<std::endl;
}

void PrintGradTable(const vector<string>& Rows,
                    const DerivMap& Derivs,
                    const SystemMap& Systems){
    
   const size_t NCols=4;
   size_t NRows=1;
   for(const auto& Di: Derivs)NRows+=(Di.second.size()/3);
   pulsar::output::Table ResultTable(NRows,NCols);
   std::array<string,NCols> ColTitles=
            {"System [Method Key]","dE/dx (a.u.)","dE/dy (a.u.)","dE/dz (a.u.)"};
   ResultTable.SetHBorder(0,'*');
   ResultTable.SetHBorder(NRows,'*');
   ResultTable.SetVBorder(1,'|');
   ResultTable.FillRow(ColTitles,0,0,NCols);
   size_t counter=0;
   for(const auto& RowI: Rows){
       ResultTable.SetHBorder(counter+1,'-');
       size_t NAtoms=Derivs.at(RowI).size()/3;
       for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
           if(AtomI==(NAtoms-NAtoms%2)/2)
               ResultTable.GetCell(counter+1,0).AddData(RowI);
           ResultTable.FillRow(&Derivs.at(RowI)[AtomI*3],++counter,1,NCols);
       }
   }
   pulsar::output::GetGlobalOut()<<ResultTable<<std::endl;
}


}//End namespace

