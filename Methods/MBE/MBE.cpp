/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <string>
#include <pulsar/math/PowerSetItr.hpp>
#include <pulsar/math/Binomial.hpp>
#include "Methods/MBE/MBE.hpp"
#include "Methods/MBE/MBECommon.hpp"
#include "Methods/MBE/MBEUtils.hpp"



using std::string;
using std::vector;
using std::map;
using std::set;
using std::stringstream;

using pulsar::modulemanager::ModuleManager;
typedef vector<double> Return_t;
typedef map<string,Return_t> DerivMap;

/*
 * TODO once computations cache results, compute interactions by
 *   "Re-running" the computations and then assembling derivatives
 */

namespace pulsarmethods{
    
    EnergyMethod::DerivReturnType MBE::Deriv_(size_t Order,const Wavefunction& Wfn){
        //Load options
        const OptionMap& DaOptions=Options();
        string MethodName=DaOptions.Get<string>("METHOD");
        string Fragmentizer=DaOptions.Get<string>("FRAGMENTIZER");

        //Make N-Mers
        const System& Mol=*Wfn.system;
        Fragmenter_t Fragger=CreateChild<SystemFragmenter>(Fragmentizer);
        SystemMap NMers=Fragger->Fragmentize(Mol);
        
        //Parse the names, bin n-mers by size, zero weights
        SNList_t SNs=BinNMers(NMers);
        map<string,double> Weights;
        for(const auto& NMerI:NMers)Weights[NMerI.first]=0.0;

        
        /* For the MBE we have the sum of all one-body interactions, the
         * sum of all two-body interactions, etc.  For each n-mer, GetCoef
         * adjusts the coefficients of the systems needed to get the energy
         * of the n-body interaction unique to that n-mer.  To get the
         * total we thus have to loop over all orders
         */
        for(size_t i=0;i<SNs.size();++i)
            for(const auto& NMerI: SNs[i])
                GetCoef(true,NMerI.first,SNs,Weights);        
        
        //Now we need the weights in the same order as the systems
        Return_t SortedWeights;
        for(const auto& NMerI: NMers)
            SortedWeights.push_back(Weights[NMerI.first]);
            
        
        //Pass info to MIM
        const string MIMName=Options().Get<string>("MIM_KEY");
        MManager().ChangeOption(MIMName,"METHODS",vector<string>({MethodName}));
        MManager().ChangeOption(MIMName,"FRAGMENTIZER",Fragmentizer);
        MManager().ChangeOption(MIMName,"WEIGHTS",SortedWeights);
        EMethod_t MIM=CreateChild<EnergyMethod>(MIMName);
        return MIM->Deriv(Order,Wfn);
    }
    

}//End namespace
