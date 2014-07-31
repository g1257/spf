
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
	typedef std::pair<RealType,RealType> PairRealRealType;

public:
	Adjustments(const EngineParamsType& engineParams)
	    : engineParams_(engineParams)
	{}

	RealType adjChemPot(const typename PsimagLite::Vector<RealType>::Type& eigs) const
	{
		typedef DensityFunctionDiag<EngineParamsType> DensityFunctionType;
		DensityFunctionType densityFunction(engineParams_,eigs);
		return adjChemPot_(densityFunction);
	}

	template<typename SomeTpemType>
	RealType adjChemPot(const typename PsimagLite::Vector<RealType>::Type& moments,
	                    const SomeTpemType& tpem) const
	{
		typedef DensityFunctionTpem<EngineParamsType,SomeTpemType> DensityFunctionType;
		DensityFunctionType densityFunction(engineParams_,moments,tpem);
		return adjChemPot_(densityFunction);
	}

private:

	template<typename SomeDensityFunctionType>
	RealType adjChemPot_(SomeDensityFunctionType& densityFunction) const
	{
		typedef PsimagLite::RootFindingBisection<SomeDensityFunctionType> RootFindingType;
		PairRealRealType aAndB = findAandB(densityFunction);
		RealType tolerance = engineParams_.adjustTolerance;
		SizeType maxIter = engineParams_.adjustMaxIter;
		RootFindingType rootFinding(densityFunction,
		                             aAndB.first,
		                             aAndB.second,
		                             maxIter,
		                             tolerance);

		RealType mu = rootFinding();
		//std::cerr<<" new mu = "<<mu<<"\n";
		return mu;
	}

	template<typename SomeDensityFunctionType>
	PairRealRealType findAandB(const SomeDensityFunctionType& function) const
	{
		for (RealType value = 1.0; value < 100.0; value++) {
			RealType value2 = function(value) * function(-value);
			if (value2 < 0) return PairRealRealType(-value,value);
		}

		throw PsimagLite::RuntimeError("RootFinding init failed\n");
	}

	const EngineParamsType& engineParams_;
}; // Adjustments
} // namespace Spf

/*@}*/
#endif// ADJUSTMENTS_H

