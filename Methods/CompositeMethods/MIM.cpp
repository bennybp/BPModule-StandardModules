#include <unordered_map>
#include <limits>
#include <algorithm>
#include <memory>

#include <pulsar/exception/Exceptions.hpp>
#include <pulsar/parallel/InitFinalize.hpp>
#include <LibTaskForce/LibTaskForce.hpp>
#include <pulsar/output/Table.hpp>
#include <pulsar/output/GlobalOutput.hpp>
#include "Methods/CompositeMethods/MIM.hpp" 

/*using pulsar::modulemanager::ModuleManager;
using LibTaskForce::Communicator;
using LibTaskForce::TaskResults;
using std::vector;
using std::map;
using std::string;

using Return_t=pulsar::modulebase::EnergyMethod::DerivReturnType;
using DerivMap=map<string,Return_t>;
typedef unsigned long ULI;*/


namespace pulsarmethods{



/*void PrintEgyTable(const vector<string>& Rows,const DerivMap& Derivs);
void PrintGradTable(const vector<string>& Rows,const DerivMap& Derivs,
                    const SystemMap& Systems);*/



pulsar::modulebase::DerivReturnType MIM::deriv_(size_t Order,const pulsar::datastore::Wavefunction& Wfn){
    return {Wfn,{0.0}};
   /*//Get the system and compute the number of degrees of freedom for the result
   const System& Mol=*Wfn.system;
   size_t DoF=1;
   for(size_t i=0;i<Order;++i)DoF*=3*Mol.size();
    
   //Establish an atom order
   AtomMap_t AtomMap=MapAtoms(Mol);
   vector<string> MethodNames=options().get<vector<string>>("METHODS");
   vector<double> Coeffs=options().get<vector<double>>("WEIGHTS");
   
   //Get the subsystems
   Fragmenter_t Fragger=create_child_from_option<SystemFragmenter>("FRAGMENTIZER");
   NMerSetType Systems=Fragger->fragmentize(Mol);   
   
   //Loop two, gettin results
   DerivMap Derivs; //Will be the derivatives per system
   Return_t TotalDeriv;
   TotalDeriv.second=vector<double>(DoF,0.0);//Final, total derivative
   
   SysI=Systems.begin();
   vector<string> RowTitles(NTasks);
   for(size_t TaskI=0;TaskI<NTasks;++TaskI){
      const System& SubSys=SysI->second.NMer;
      std::stringstream ss;
      for(const std::string& Name: SysI->first)
        ss<<Name<<" ";
      ss<<"["<<MethodNames[SameMethod?0:TaskI]<<"]";
      Derivs[ss.str()]=Results[TaskI];
      RowTitles[TaskI]=ss.str();
      Fillderiv(TotalDeriv.second,Results[TaskI].second,
                SysI->second.Weight*Coeffs[TaskI],SubSys,
                AtomMap,MapAtoms(SubSys),Order);
      if(!SameSystem)++SysI;
   }
   if(Order==0)PrintEgyTable(RowTitles,Derivs);
   else if(Order==1)PrintGradTable(RowTitles,Derivs);
   return TotalDeriv;*/
}

/*void PrintEgyTable(const vector<string>& Rows,const DerivMap& Derivs){
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
         ResultTable.FillRow(Derivs.at(Rows[i]).second.data(),i+1,1,NCols);
   pulsar::output::GetGlobalOut()<<ResultTable<<std::endl;
}

void PrintGradTable(const vector<string>& Rows,
                    const DerivMap& Derivs){
    
   const size_t NCols=4;
   size_t NRows=1;
   for(const auto& Di: Derivs)NRows+=(Di.second.second.size()/3);
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
       size_t NAtoms=Derivs.at(RowI).second.size()/3;
       for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
           if(AtomI==(NAtoms-NAtoms%2)/2)
               ResultTable.GetCell(counter+1,0).AddData(RowI);
           ResultTable.FillRow(&Derivs.at(RowI).second[AtomI*3],++counter,1,NCols);
       }
   }
   pulsar::output::GetGlobalOut()<<ResultTable<<std::endl;
}*/


}//End namespace

