
/** \ingroup SPF */
/*@{*/

/*! \file Adjustments.h
 *
 *  Adjustments such as chemical potential
 *
 */
#ifndef ADJUSTMENTS_H
#define ADJUSTMENTS_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite
#include "DensityFunction.h"
#include "RootFindingNewton.h"
#include "RootFindingBisection.h"

namespace Spf {
	template<typename EngineParamsType>
	class Adjustments {
		typedef typename EngineParamsType::RealType RealType;
		typedef DensityFunction<EngineParamsType> DensityFunctionType;
		// typedef RootFindingNewton<DensityFunctionType> RootFindingType;
		typedef RootFindingBisection<DensityFunctionType> RootFindingType;
	public:
		Adjustments(const EngineParamsType& engineParams,
		            size_t maxIter=100,
		            RealType tolerance=1.0e-3)
		: engineParams_(engineParams),maxIter_(maxIter),tolerance_(tolerance)
		{
		}
		
		RealType adjChemPot(const std::vector<RealType>& eigs) const
		{
			DensityFunctionType densityFunction(engineParams_,eigs);
			RootFindingType  rootFinding(densityFunction);
			
			RealType mu=engineParams_.mu;
			//std::cerr<<"Old mu="<<mu<<" ";
			rootFinding(mu);
			//std::cerr<<" new mu = "<<mu<<"\n";
			return mu;
		}
		
		/* template<typename SomeOutputType>
		void print(SomeOutputType& os) const
		{
			os<<"Adjustments: mu="<<engineParams_.mu<<"\n";
		}*/
		
	private:
		
		const EngineParamsType& engineParams_;
		size_t maxIter_;
		RealType tolerance_;
	}; // Adjustments
} // namespace Spf

/*@}*/
#endif// ADJUSTMENTS_H

