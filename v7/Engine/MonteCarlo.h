
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
		typedef std::pair<size_t,size_t> PairType;
		
	public:
		typedef typename EngineParamsType::FieldType FieldType;
		
		MonteCarlo(const EngineParamsType& engineParams,ModelType& model,AlgorithmType& algorithm) 
			: engineParams_(engineParams),model_(model),rng_(),algorithm_(algorithm) { }
		
		template<typename DynVarsType> // DynVarsType is either Spin or Phonon
		PairType operator()(DynVarsType& dynVars, size_t iter)
		{
			PairType acc = PairType(0,0);
			model_.set(dynVars);
			algorithm_.init();
			for (size_t j=0;j<dynVars.size;j++) {
				size_t i = model_.proposeSite(j,rng_);	
				model_.proposeChange(i,rng_);
				
				bool flag= algorithm_.isAccepted(i);
				//std::cerr<<"flag="<<flag<<"\n";
				if (flag && !dynVars.isFrozen) { // Accepted
					algorithm_.accept(i);
					acc.first++;
				} else { // not accepted
					//engineParams_.mu=oldmu;
				}
				acc.second++;
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

