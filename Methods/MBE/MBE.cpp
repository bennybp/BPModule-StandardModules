/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <limits>
#include <bpmodule/system/System.hpp>
#include <bpmodule/modulebase/SystemFragmenter.hpp>
#include <bpmodule/datastore/OptionMap.hpp>
#include <bpmodule/parallel/InitFinalize.hpp>
#include <bpmodule/math/PowerSetItr.hpp>
#include <bpmodule/math/Binomial.hpp>
#include <bpmodule/output/Table.hpp>
#include <LibTaskForce.hpp>
#include "Methods/MBE/MBE.hpp"

using bpmodule::datastore::OptionMap;
using bpmodule::system::Atom;
using bpmodule::system::System;
using bpmodule::modulemanager::ModuleManager;
using bpmodule::modulebase::SystemFragmenter;
using bpmodule::modulebase::EnergyMethod;
using bpmodule::system::SystemMap;
using LibTaskForce::Communicator;
using LibTaskForce::TaskResults;


typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;
typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;
typedef std::vector<double> Return_t;
typedef std::vector<std::string> String_t;
typedef std::map<std::string,Return_t> DerivMap;

/*
 * TODO when I can compile and know I didn't break things:
 * - Remove task iteration and derivative putting together
 *   and write in terms of MIM
 * - Make N-Mer formation an option of the system fragmenter
 *   - Put it in base fragmenter
 * - Once computations cache results, compute interactions by
 *   "Re-running" the computations and then assembling derivatives
 */

namespace bpmethods{

    class Task{
    private:
        ModuleManager& MM_;
        std::string Key_;
        unsigned long ID_;
        const System& Sys_;  //! \todo Could be made a shared pointer,
                             //        so it isn't copied in operator()?
                             // ie, Constructor would have to take a shared pointer
    public:
        Task(ModuleManager& MM,
             std::string Key,
             unsigned long ID,
             const System& Sys):
                MM_(MM),Key_(Key),ID_(ID),Sys_(Sys){ }
        Return_t operator()(size_t DerivOrder)const{
            EMethod_t Method=MM_.GetModule<EnergyMethod>(Key_,ID_);
            Method->Wfn().system = std::make_shared<const System>(Sys_);
            return Method->Deriv(DerivOrder);
        }
    };

    //Computes n-body interaction by recursion
    Return_t NBodyEgy(size_t DoF,
                      const std::vector<DerivMap>& Ints,
                      const std::vector<DerivMap>& Derivs,
                      const String_t& Frags);
    String_t split(const std::string &s){
        String_t elems;
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss,item,'_'))elems.push_back(item);
        return elems;
    }
    Return_t MBE::DerivImpl(size_t Order)const{
        //Load options
        const OptionMap& DaOptions=Options();
        std::string MethodName=DaOptions.Get<std::string>("METHOD");
        size_t N=DaOptions.Get<size_t>("TRUNCATION_ORDER");
        std::map<size_t,double> Truncs=
                DaOptions.Get<std::map<size_t,double>>("DISTANCE_THRESHOLDS");

        //Compute size of super system basis
        const System& Mol=*(Wfn().system);
        size_t DoF=1;
        for(size_t i=0;i<Order;++i)DoF*=3*Mol.Size();


        //Make N-Mers
        Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>("FRAG");
        std::vector<SystemMap> NMers(N);
        NMers[0]=Fragger->Fragmentize(Mol);
        for(size_t n=2;n<=N;++n){
            //Defaults to about 1e308 a.u., or about 1e271 times the size of
            //the observable universe...I think that is equivalent to no cut-off
            double dist=(Truncs.count(n)==1?
                         Truncs.at(n):std::numeric_limits<double>::max());
            NMers[n-1]=bpmodule::system::MakeNMers(NMers[0],n,dist);
        }

        //Run N-Mers
        const Communicator& Comm=bpmodule::parallel::GetEnv().Comm();
        Communicator NewComm=Comm.Split(Comm.NThreads(),1);
        TaskResults<Return_t> Results(NewComm);
        for(size_t n=N;n>0;--n){
            for(const auto& NMer:NMers[n-1]){
                Results.push_back(NewComm.AddTask(
                                                  Task(MManager(),MethodName,ID(),NMer.second),Order));
            }
        }

        //Get results (we loop in same order)
        std::vector<DerivMap> Derivs(N);//Each individual result
        size_t counter=0;
        for(size_t n=N;n>0;--n)
            for(const auto& NMer:NMers[n-1]){
                Derivs[n-1].emplace(NMer.first,Results[counter++]);
                //\todo Expand into supersystem basis
            }

        //Compute interactions (and derivatives of interactions)
        //All quantities need to be in supersystem basis for this to work
        std::vector<Return_t> TotalEgysByOrder(N,Return_t(DoF,0.0));
        std::vector<DerivMap> Interactions(N);
        for(size_t n=0;n<N;++n){
            for(const auto& NMer:NMers[n]){
                for(size_t i=0;i<DoF;++i)
                    TotalEgysByOrder[n][i]+=Derivs[n][NMer.first][i];
                Interactions[n].emplace(NMer.first,
                                        NBodyEgy(DoF,Interactions,Derivs,split(NMer.first))
                                        );
            }
        }

        //Total by summing up all interactions
        std::vector<Return_t> Egy(N,Return_t(DoF,0.0));
        std::vector<Return_t> Egy2(N,Return_t(DoF,0.0));

        /* We total two ways:
         *  
         * 1. Sum of all interactions
         * 2. Closed form
         *
         * The reason is the closed form should propagate error less, but
         * in general we will need sum of interactions.  This provides us
         * a manner to establish accuracy.
         * 
         * Note on 2: if n==N Binomial coefficients push implementations b/c
         * it involves negative binomial coefficients and ones for which the
         * "choose k" is greater.  The expansion is easy in that case though,
         * it's just the supersystem energy.
         */
        const size_t NFrags=NMers[0].size();
        for(size_t BigOrder=1;BigOrder<=N;++BigOrder){
            if(BigOrder>1)Egy2[BigOrder-1]=Egy2[BigOrder-2];
            for(const auto& Int:Interactions[BigOrder-1])
                for(size_t i=0;i<DoF;++i)
                    Egy2[BigOrder-1][i]+=Int.second[i];
            if(BigOrder<N){
                for(size_t n=0;n<BigOrder;++n){
                    double coef=pow(-1,BigOrder-n-1)*
                            bpmodule::math::BinomialCoefficient(NFrags-n-2,BigOrder-n-1);
                    for(size_t i=0;i<DoF;++i)
                        Egy[BigOrder-1][i]+=coef*TotalEgysByOrder[n][i];
                }
            }else Egy[BigOrder-1]=TotalEgysByOrder[BigOrder-1];
        }
        /*for(size_t n=0;n<N;n++){
            for(size_t i=0;i<DoF;++i)
                std::cout<<Egy[n][i]<<" "<<Egy2[n][i]<<std::endl;
        }*/
        bpmodule::output::Table TableOfResults(N+1,3);
        std::array<std::string,3> ColTitles={"","Closed-Form","SumOfInts"};
        std::vector<std::string> RowTitles(N);
        Return_t Col1,Col2;
        for(size_t i=0;i<N;++i){
            std::stringstream ss;
            ss<<i+1<<"-body energy(a.u.):";
            RowTitles[i]=ss.str();
            Col1.push_back(Egy[i][0]);
            Col2.push_back(Egy2[i][0]);
        }
        
        TableOfResults.FillRow(ColTitles,0,0,3);
        TableOfResults.FillCol(RowTitles,0,1,N+1);
        TableOfResults.FillCol(Col1,1,1,N+1);
        TableOfResults.FillCol(Col2,2,1,N+1);
        TableOfResults.SetHBorder(0,'*');
        TableOfResults.SetHBorder(1,'-');
        TableOfResults.SetHBorder(N+1,'*');
        TableOfResults.SetVBorder(0,'|');
        TableOfResults.SetVBorder(1,'|');
        TableOfResults.SetVBorder(3,'|');
        TableOfResults.SetAlign(0,bpmodule::output::LEFT);
        std::cout<<TableOfResults<<std::endl;
        return Egy[N-1];
    }
    /* This function computes the n-body interaction of a given n-mer, recursively.
     * Say we want the four-body interaction among IJKL, that's the energy
     * of IJKL less the interactions: IJK, IJL, JKL, IJ, IK, IL, JK, JL, KL, I, J,
     * K, L, which is the power set of IJKL w/o the empty set
     */
    Return_t NBodyEgy(size_t DoF,
                      const std::vector<DerivMap>& Ints,
                      const std::vector<DerivMap>& Derivs,
                      const String_t& Frags){
        std::stringstream ss;
        ss<<Frags[0];
        const size_t n=Frags.size();
        for(size_t i=1;i<n;++i)ss<<"_"<<Frags[i];

        //If we already computed this interaction just return
        if(Ints[n-1].count(ss.str())==1)return Ints[n-1].at(ss.str());

        //Didn't compute the energy? Assume the interaction is negligible
        if(Derivs[n-1].count(ss.str())==0)return Return_t();

        //Otherwise, we're computing it
        Return_t Egy=Derivs[n-1].at(ss.str());

        bpmodule::math::PowerSetItr<String_t> SubInts(Frags,1,n-1);
        while(!SubInts.Done()){
            Return_t TempEgy=NBodyEgy(DoF,Ints,Derivs,*SubInts);
            if(TempEgy.size()>0)//Don't subtract out 0 interactions
                for(size_t i=0;i<DoF;++i)Egy[i]-=TempEgy[i];
            ++SubInts;
        }
        return Egy;
    }

}//End namespace
