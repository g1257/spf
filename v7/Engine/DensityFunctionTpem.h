
/** \ingroup SPF */
/*@{*/

/*! \file DensityFunctionTpem.h
 *
 *  DensityFunctionTpem such as chemical potential
 *
 */
#ifndef DENSITY_FUNCTION_TPEM_H
#define DENSITY_FUNCTION_TPEM_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/

namespace Spf {
	template<typename EngineParamsType,typename TpemType>
	class DensityFunctionTpem {
	public:
		typedef typename EngineParamsType::RealType RealType;
		typedef typename TpemType::NumberFunctorType NumberFunctorType;
		typedef typename TpemType::TpemParametersType TpemParamtersType;

		DensityFunctionTpem(const EngineParamsType& engineParams,
		                    const std::vector<RealType>& moments,
		                    const TpemType& tpem)
		: engineParams_(engineParams),
		  tpem_(tpem),
		  moments_(moments),
		  numberFunctor_(tpem_.tpemParameters())
		{}

		RealType operator()(const RealType& mu) const
		{
			RealType muSaved = engineParams_.mu;
			engineParams_.mu = mu;
			std::vector<RealType> numberCoeffs(moments_.size());
			tpem_.calcCoeffs(numberCoeffs,numberFunctor_);
			engineParams_.mu = muSaved;
			return tpem_.expand(moments_, numberCoeffs) - engineParams_.carriers;
		}

		// Derivative of n(mu) with respect to mu
		RealType derivative(const RealType& mu) const
		{
			throw std::runtime_error("Unimplemented yet\n");
		}
		
		const EngineParamsType& engineParams_;
		const TpemType& tpem_;
		const std::vector<RealType>& moments_;
		NumberFunctorType numberFunctor_;
	}; // DensityFunctionTpem
} // namespace Spf

/*@}*/
#endif// DENSITY_FUNCTION_TPEM_H

