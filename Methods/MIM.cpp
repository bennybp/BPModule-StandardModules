#include "Methods/MIM/MIM.hpp" //Note that MIM.hpp is not in a folder yet b/c I can't figure out how to make folders on GitHub

using bpmodule::modulemanager::ModuleManager;
using bpmodule::system::System;

typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;
typedef std::vector<double> Return_t;

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


Return_t MIM::DerivImpl(size_t Order)const{
   const OptionMap& DaOptions=Options();
   std::vector<std::string> MethodNames=DaOptions.Get<std::vector<std::string>>("METHODS");
   std::vector<bpmodule::system::System> Systems;//TODO: figure out how to get these, can't pass through options
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
   for(size_t TaskI=0; TaskI<NTasks; ++TaskI){
      Results.push_back(
         NewComm.AddTask(
            Task(MManager(),
                 SameMethods?MethodNames[0]:MethodNames[TaskI],
                 ID(),
                 SameSystem?Systems[0]:Systems[TaskI];
            ),
         Order);
   }
   
   
   return Return_t(0.0);
}

}//End namespace

