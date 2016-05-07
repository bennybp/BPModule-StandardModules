/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <sstream>
#include <pulsar/math/PowerSetItr.hpp>
#include "Methods/BSSE/BSSEKernel.hpp"
#include "Methods/MBE/MBECommon.hpp"

namespace pulsarmethods{
    using std::string;
    using std::vector;
    using std::map;
    
    typedef vector<double> Return_t;

   inline System SwitchBasis(const System& OldRealSys,const System& NewU){
      return NewU.Partition([&OldRealSys](const Atom& AtomI){
                            return OldRealSys.Contains(AtomI);});
   }
    
RealGhostData GhostTheSystem(const pulsar::system::System& Sys){
    RealGhostData Data;
    
   //We use two loops so that the real atoms come first and then the ghosts
   //ultimately this will make derivatives easier...
   for(const Atom& AtomI: Sys){
       Data.NewSystem<<AtomI;
       Data.Atom2Idx.emplace(AtomI,Data.Atom2Idx.size());
   }
   size_t NAtoms=Sys.Size();
   Data.Atom2RealIdx=Data.Atom2Idx;
   for(const Atom& AtomI: Sys){
       Atom Ghost=MakeGhost(NAtoms++,AtomI); 
       Data.Real2Ghost.emplace(AtomI,Ghost);
       Data.NewSystem<<Ghost;
       Data.Atom2Idx.emplace(Ghost,Data.Atom2Idx.size());
       Data.Atom2RealIdx.emplace(Ghost,Data.Atom2Idx.at(AtomI));
    }
   Data.RealSystem=std::unique_ptr<System>(new System(Data.NewSystem,false));
   for(const Atom& AtomI: Sys)Data.RealSystem->Insert(AtomI);
   return Data;
}
   inline SN_t SetDiff(const SN_t& FullSN,const SN_t& SN){
        SN_t Ghosts;
        std::set_difference(FullSN.begin(),FullSN.end(),
                            SN.begin(),SN.end(),
                            std::inserter(Ghosts,Ghosts.begin()));
        return Ghosts;
   }  

   string MakeName(const SN_t& FullSN,const SN_t& SN){
        SN_t NewSN=SN,Ghosts=SetDiff(FullSN,SN);
        std::stringstream ss;
        for(size_t i: Ghosts)ss<<"Gh("<<i<<")_";
        size_t counter=0;
        for(size_t i: SN){
            ss<<i;
            if(++counter<SN.size())ss<<"_";
        }
        return ss.str();
    }
   
   std::map<string,double> SSFCKernel(SystemMap& NMers,
                                      const RealGhostData& NewSys, 
                                      size_t MinRealOrder, 
                                      size_t MaxRealOrder,
                                      size_t MinGhostOrder,
                                      size_t MaxGhostOrder){
        /*  Note that the GetCoef function from
         *  MBEUtil will not weight subsystems that have no serial number in the
         *  function's input.  Hence we don't add ghost only systems to AllSNs
         *  and we don't have to worry about them manifesting.
         */
        SNList_t AllSNs=BinNMers(NMers),RealSNs=BinNMers(NMers);
        std::map<string,double> Coeffs;
        System NewU(NewSys.NewSystem,true);
        System Empty=NewU.Partition([](const Atom&){return false;});
        System RealSystem=SwitchBasis(*NewSys.RealSystem,NewU);
        //Get the serial number of the full system and make ghost monomers
        SN_t FullSN;
        size_t NFrags=RealSNs[0].size();
        std::map<size_t,System> GhostMonos;
        for(size_t i=0;i<NFrags;++i){
            FullSN.insert(i);
            GhostMonos.emplace(i,Empty);
            for(const Atom& AtomI:NMers.at(RealSNs[0][{i}]))
                GhostMonos.at(i).Insert(NewSys.Real2Ghost.at(AtomI));
        }
        
        ///By going in order we ensure smallest systems are added first thus
        ///when recursion occurs they will be there
        for(size_t i=MinRealOrder-1;i<MaxRealOrder;++i){
            for(const auto& NMerI: RealSNs[i]){//Loop over real nmers
                const string& OldName=NMerI.second;
                System NewFrag=SwitchBasis(NMers.at(OldName),NewU);
                //Get only frags that can be ghosts
                SN_t EffectiveFullSN=SetDiff(FullSN,NMerI.first);
                pulsar::math::PowerSetItr<SN_t> 
                    Basis(EffectiveFullSN,MinGhostOrder,
                          std::min(MaxGhostOrder,EffectiveFullSN.size()));
                
                //For each n-mer need to add appropriate combinations of ghosts
                while(!Basis.Done()){
                    SN_t SN=NMerI.first;
                    SN_t Ghosts=*Basis;
                    for(size_t Gi:Ghosts){
                        NewFrag+=GhostMonos.at(Gi);
                       //We offset ghost frags by NFrags to ensure they don't
                       //get in the way
                        SN.insert(Gi+NFrags);
                    }
                    string NewName=MakeName(SN,SN);
                    if(AllSNs.size()<SN.size())AllSNs.resize(SN.size());
                    AllSNs[SN.size()-1][SN]=NewName;
                    GetCoef(false,SN,AllSNs,Coeffs);
                    NMers.emplace(NewName,NewFrag);
                    ++Basis;
                }
            }
        }
        return Coeffs;
   }
   
    
   Return_t RunCalcs(const pulsar::system::SystemMap& AllFrags,
                           const map<string,double>& Coeffs,
                           const RealGhostData& Data,
                           size_t Order,
                           ID_t ID,
                           pulsar::modulemanager::ModuleManager& MM,
                           const string& MethodName,
                           const string& MIMName){
        size_t NAtoms=Data.RealSystem->Size();
        size_t DoF=1;
        for(size_t i=0;i<Order;++i)DoF*=3*NAtoms;
       
        System NewU(Data.NewSystem,true); 
        Return_t Cs;
        vector<string> Names;
        vector<size_t> FragAtoms,Offsets;
        for(const auto& NMerI: AllFrags){
            const string& Name=NMerI.first;
            const System& Frag=NMerI.second;
            Names.push_back(Name);
            Cs.push_back(Coeffs.at(Name));
            Offsets.push_back(Frag.Size());
            for(const Atom& AtomI: Frag)
                FragAtoms.push_back(Data.Atom2Idx.at(AtomI));
        }
        
        //MM.ReplaceKey("BSSE_FRAG","UserDefined");
        string MIMKey=MM.GenerateUniqueKey();
        string FragKey=MM.GenerateUniqueKey();
        MM.DuplicateKey(MIMName,MIMKey);
        MM.DuplicateKey("BP_UD_FRAG",FragKey);
        //TODO: get rid of dependence on BP_UD_FRAG
        MM.ChangeOption(FragKey,"FRAGMENT_NAMES",Names);
        MM.ChangeOption(FragKey,"ATOMS_PER_FRAG",Offsets);
        MM.ChangeOption(FragKey,"FRAGMENTS",FragAtoms);
        MM.ChangeOption(MIMKey,"FRAGMENTIZER",FragKey);
        MM.ChangeOption(MIMKey,"METHODS",vector<string>({MethodName}));
        MM.ChangeOption(MIMKey,"WEIGHTS",Cs);

        
        EMethod_t MIM=MM.GetModule<EnergyMethod>(MIMKey,ID);
        MIM->InitialWfn().GetSystem()=std::make_shared<System>(Data.NewSystem,true);
        Return_t FullDeriv=MIM->Deriv(Order);
        Return_t FinalDeriv(DoF,0.0);
        FillDeriv(FinalDeriv,FullDeriv,1.0,NewU,Data.Atom2RealIdx,
                                                Data.Atom2Idx,Order);
        return FinalDeriv;
   }
   
}//End namespace
