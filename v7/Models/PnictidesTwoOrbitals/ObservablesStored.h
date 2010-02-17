
/** \ingroup SPF */
/*@{*/

/*! \file ObservablesStored.h
 *
 *  ObservablesStored for pnictidesTwoOrbitals model
 *
 */

#ifndef OBSERVABLES_STORED_H
#define OBSERVABLES_STORED_H
#include "Vector.h"

namespace Spf {
	template<typename SpinOperationsType>
	class ObservablesStored {
		
		typedef typename SpinOperationsType::DynVarsType DynVarsType;
		typedef typename DynVarsType::FieldType FieldType;
		typedef PsimagLite::Vector<FieldType> VectorType;
		
	public:
		ObservablesStored(SpinOperationsType& spinOperations,size_t vol,size_t nbands) : 
			spinOperations_(spinOperations),cc_(vol),lc_(nbands*vol)
		{}
				
		template<typename GreenFunctionType>
		void operator()(const DynVarsType& spins,GreenFunctionType& greenFunction)
		{
			spinOperations_.classicalCorrelations(cc_,spins);
			greenFunction.localCharge(lc_);
			counter_++;
		}
		
		void finalize(std::ostream& fout)
		{
			cc_/=counter_;
			fout<<"#ClassicalCorrelations:\n";
			fout<<cc_;	
		}
		
	private:
		SpinOperationsType& spinOperations_;
		VectorType cc_;
		VectorType lc_;
		size_t counter_;
	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBSERVABLES_STORED_H
