#ifndef MIM_HPP_
#define MIM_HPP_

/* Disclaimer:  Ryan wrote this on GitHub, while trying to kill time.  It likely does not compile*/

#include <vector>
#include <bpmodule/modulebase/EnergyMethod.hpp>

namespace bpmethods{

/** \brief This is a class that implements the Molecules in Molecules (MIM) method of Mayhall and Rahavachari (sp)
 *
 *  For all intents and purposes, you can think of MIM as ONIOM where each computation can also be a MBE.  Why
 *  that required a new name is beyond me, but given that I couldn't think of a better name for this class I
 *  accepted that.  Basically this class implements all composite methods that can be written in the form:
 *  \f[
 *     E=\sum_{i=1}^{M}c_iE_i,
 *  \f]
 *  where \f$E\f$ is the total, approximate energy which is obtained as a linear combination of \f$M\f$ other
 *  energies.  For example, the standard trick to compute the CCSD(T) energy in a large basis set:
 *  \f[
 *     E_{CCSD(T)}^{Large}=E_{MP2}^{Large}+\left(E_{CCSD(T)}^{Small}-E_{MP2}^{Small}\right)
 *  \f]
 *  is a maifestation of this, as are QM/MM methods, the MBE, and BSSE corrections.
 *
 *  In particular, this class factors out the code required to submit the various computations in parallel and will
 *  then also handle "stitching" the results together.
 *
 *  This class needs three things from you: a list of systems to run, the weight of each system in the resulting expansion,
 *  and the keys for each method that should be run.  Respectively, these options are SYSTEMS, WEIGHTS, and METHODS.
 */
class MIM : public bpmodule::modulebase::EnergyMethod {
   private:
      ///The type of the base class
      typedef bpmodule::modulebase::EnergyMethod Base_t;
      ///The function that actually implements the derivative computation
      std::vector<double> DerivImpl(size_t Order)const;
   public:
      ///Pulls in the base class's methods
      using Base_t::EnergyMethod;
      
      ///The method that the base class will actually call
      std::vector<double> Deriv_(size_t Order){
         if(Order>=1)return Base_t::Deriv_(Order);
         return DerivImpl(Order);
      }
};

}//End namespace

#endif /*End header guard MIM_HPP_ */
