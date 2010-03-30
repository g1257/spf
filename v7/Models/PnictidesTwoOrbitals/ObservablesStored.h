
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
		ObservablesStored(SpinOperationsType& spinOperations,size_t vol,size_t dof) : 
			spinOperations_(spinOperations),cc_(vol),lc_(dof*vol),counter_(0)
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
			divideAndPrint(fout,cc_,"#ClassicalCorrelations:");
			divideAndPrint(fout,lc_,"#LocalCharge:");
		}
		
	private:
		SpinOperationsType& spinOperations_;
		VectorType cc_;
		VectorType lc_;
		size_t counter_;
		
		void divideAndPrint(std::ostream& fout,VectorType& v,const std::string& label)
		{
			v /= counter_;
			fout<<label<<"\n";
			fout<<v;
		}

	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBSERVABLES_STORED_H
