
/** \ingroup SPF */
/*@{*/

/*! \file DensityFunction.h
 *
 *  DensityFunction such as chemical potential
 *
 */
#ifndef DENSITY_FUNCTION_H
#define DENSITY_FUNCTION_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/

namespace Spf {
	template<typename EngineParamsType>
	class DensityFunction {
	public:
		typedef typename EngineParamsType::FieldType FieldType;
		
		DensityFunction(const EngineParamsType& engineParams,
						const std::vector<FieldType>& eig)
		: engineParams_(engineParams),eig_(eig)
		{
		}
		
		FieldType operator()(const FieldType& mu) const
		{
			FieldType sum = 0;
			for (size_t i=0;i<eig_.size();i++)
				sum += PsimagLite::fermi((eig_[i]-mu)*engineParams_.beta);
			return sum - engineParams_.carriers;
		}

		// Derivative of n(mu) with respect to mu
		FieldType derivative(const FieldType& mu) const
		{
			FieldType sum=0;
			const FieldType& beta = engineParams_.beta;
			for (size_t i=0;i<eig_.size();i++)
				sum -= PsimagLite::fermiPrime((eig_[i]-mu)*beta)*beta;
			return sum;
		}
		
		const EngineParamsType& engineParams_;
		const std::vector<FieldType>& eig_;
	}; // DensityFunction
} // namespace Spf

/*@}*/
#endif// DENSITY_FUNCTION_H

