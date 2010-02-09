
/** \ingroup SPF */
/*@{*/

/*! \file MonteCarlo.h
 *
 *  Monte Carlo for SPF
 *
 */
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include "Utils.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename EngineParamsType,typename ModelType,typename AlgorithmType,typename RandomNumberGeneratorType>
	class MonteCarlo {
		typedef typename EngineParamsType::FieldType FieldType;
		
	public:
		MonteCarlo(const EngineParamsType& engineParams,ModelType& model,AlgorithmType& algorithm) 
			: engineParams_(engineParams),model_(model),rng_(),algorithm_(algorithm) { }
		
		template<typename DynVarsType>
		size_t operator()(DynVarsType& dynVars, size_t iter)
		{
			size_t acc = 0;
			algorithm_.init(dynVars);
			for (size_t i=0;i<dynVars.size();i++) {
					
				model_.propose(i,rng_);
				
				bool flag= algorithm_.isAccepted(i);
					
				if (flag && !dynVars.isFrozen) { // Accepted
					algorithm_.accept(i);
					acc++;
				} else { // not accepted
					//engineParams_.mu=oldmu;
				}
			} // lattice sweep
			return acc;
		}
	private:
		
		
		const EngineParamsType& engineParams_;
		ModelType& model_;
		RandomNumberGeneratorType rng_;
		AlgorithmType& algorithm_;
	}; // MonteCarlo
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H

