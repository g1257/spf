
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
#include "DensityFunctionDiag.h"
#include "DensityFunctionTpem.h"
#include "RootFindingNewton.h"
#include "RootFindingBisection.h"

namespace Spf {
	template<typename EngineParamsType>
	class Adjustments {
		typedef typename EngineParamsType::RealType RealType;
	public:
		Adjustments(const EngineParamsType& engineParams)
		: engineParams_(engineParams)
		{}

		RealType adjChemPot(const typename PsimagLite::Vector<RealType>::Type& eigs) const
		{
			typedef DensityFunctionDiag<EngineParamsType> DensityFunctionType;
			DensityFunctionType densityFunction(engineParams_,eigs);
			typedef PsimagLite::RootFindingBisection<DensityFunctionType> RootFindingType;
			// typedef RootFindingNewton<DensityFunctionType> RootFindingType;
			RealType tolerance = engineParams_.adjustTolerance;
			SizeType maxIter = engineParams_.adjustMaxIter;
			RootFindingType  rootFinding(densityFunction,-10,10,maxIter,tolerance);

			RealType mu = rootFinding();
			//std::cerr<<" new mu = "<<mu<<"\n";
			return mu;
		}

		template<typename SomeTpemType>
		RealType adjChemPot(const typename PsimagLite::Vector<RealType>::Type& moments,
		                    const SomeTpemType& tpem) const
		{
			typedef DensityFunctionTpem<EngineParamsType,SomeTpemType> DensityFunctionType;
			DensityFunctionType densityFunction(engineParams_,moments,tpem);
			typedef PsimagLite::RootFindingBisection<DensityFunctionType> RootFindingType;
			// typedef RootFindingNewton<DensityFunctionType> RootFindingType;
			RealType tolerance = engineParams_.adjustTolerance;
			SizeType maxIter = engineParams_.adjustMaxIter;
			RootFindingType  rootFinding(densityFunction,-10,10,maxIter,tolerance);

			RealType mu = rootFinding();
			//std::cerr<<" new mu = "<<mu<<"\n";
			return mu;
		}

	private:

		const EngineParamsType& engineParams_;
	}; // Adjustments
} // namespace Spf

/*@}*/
#endif// ADJUSTMENTS_H

