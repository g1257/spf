
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
		ObservablesStored(SpinOperationsType& spinOperations,size_t vol) : 
			spinOperations_(spinOperations),cc_(vol)
		{}
				
		void operator()(const DynVarsType& spins)
		{
			spinOperations_.classicalCorrelations(cc_,spins);
			counter_++;
		}
		
		void finalize(std::ostream& fout)
		{
			//cc_/=counter_;
			fout<<"#ClassicalCorrelations:\n";
			fout<<cc_;	
		}
		
	private:
		SpinOperationsType& spinOperations_;
		VectorType cc_;
		size_t counter_;
	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBSERVABLES_STORED_H
