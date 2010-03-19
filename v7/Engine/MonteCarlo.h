
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
	template<typename EngineParamsType,typename OperationsType,typename AlgorithmType,typename RandomNumberGeneratorType,
 typename DynVarsType>
	class MonteCarlo {
		typedef std::pair<size_t,size_t> PairType;
		
	public:
		typedef typename EngineParamsType::FieldType FieldType;
		
		MonteCarlo(const EngineParamsType& engineParams,OperationsType& ops,AlgorithmType& algorithm,
			   RandomNumberGeneratorType& rng) 
			: engineParams_(engineParams),ops_(ops),rng_(rng),algorithm_(algorithm) { }
		
		PairType operator()(DynVarsType& dynVars, size_t iter)
		{
			PairType acc = PairType(0,0);
			ops_.set(dynVars);
			algorithm_.init();
			//std::cerr<<"F:"<<rng_()<<"\n";
			for (size_t j=0;j<dynVars.size;j++) {
				size_t i = ops_.proposeSite(j,rng_);	
				ops_.proposeChange(i,rng_);
				//FieldType oldmu = engineParams_.mu;
				bool flag= algorithm_.isAccepted(i);
				//std::cerr<<"New mu="<<engineParams_.mu<<"\n";
				//std::cerr<<"flag="<<flag<<"\n";
				if (flag && !dynVars.isFrozen) { // Accepted
					algorithm_.accept(i);
					acc.first++;
				} else { // not accepted
					//engineParams_.mu=oldmu;
					//std::cerr<<"Not accepted: oldmu="<<oldmu<<"\n";
				}
				acc.second++;
			} // lattice sweep
			return acc;
		}
		
	private:
		
		const EngineParamsType& engineParams_;
		OperationsType& ops_;
		RandomNumberGeneratorType& rng_;
		AlgorithmType& algorithm_;
		
	}; // MonteCarlo
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H

