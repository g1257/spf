
/** \ingroup SPF */
/*@{*/

/*! \file DensityFunctionDiag.h
 *
 *  DensityFunctionDiag such as chemical potential
 *
 */
#ifndef DENSITY_FUNCTION_DIAG_H
#define DENSITY_FUNCTION_DIAG_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/

namespace Spf {
	template<typename EngineParamsType>
	class DensityFunctionDiag {
	public:
		typedef typename EngineParamsType::RealType RealType;
		
		DensityFunctionDiag(const EngineParamsType& engineParams,
						const std::vector<RealType>& eig)
		: engineParams_(engineParams),eig_(eig)
		{
		}

		RealType operator()(const RealType& mu) const
		{
			RealType sum = 0;
			for (size_t i=0;i<eig_.size();i++)
				sum += PsimagLite::fermi((eig_[i]-mu)*engineParams_.beta);
			return sum - engineParams_.carriers;
		}

		// Derivative of n(mu) with respect to mu
		RealType derivative(const RealType& mu) const
		{
			RealType sum=0;
			const RealType& beta = engineParams_.beta;
			for (size_t i=0;i<eig_.size();i++)
				sum -= PsimagLite::fermiPrime((eig_[i]-mu)*beta)*beta;
			return sum;
		}
		
		const EngineParamsType& engineParams_;
		const std::vector<RealType>& eig_;
	}; // DensityFunctionDiag
} // namespace Spf

/*@}*/
#endif// DENSITY_FUNCTION_DIAG_H

