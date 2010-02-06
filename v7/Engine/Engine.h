
/** \ingroup SPF */
/*@{*/

/*! \file Engine.h
 *
 *  SPF Engine
 *
 */
#ifndef SPF_ENGINE_H
#define SPF_ENGINE_H
#include "Utils.h"

namespace Spf {
	
	template<typename ParametersType,typename ModelType,typename ConcurrencyType>
	class Engine {
		
		typedef typename ModelType::DynVarsType DynVarsType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,DynVarsType& dynVars,ConcurrencyType& concurrency) 
			: params_(params),model_(model),dynVars_(dynVars),concurrency_(concurrency)
		{
		}
				
		void main()
		{
			thermalize();
			// announce thermalization done
			measure();
			// announce measurements done
			//finalize(); FIXME
		}
		
		private:
		
		void thermalize()
		{
			for (size_t iter=0;iter<params_.iterTherm;iter++) {
				utils::printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				model_.doMonteCarlo(dynVars_,iter);
			}
		}
		
		void measure()
		{
			for (size_t iter=0;iter<params_.iterEffective;iter++) {
				utils::printProgress(iter,params_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					model_.doMonteCarlo(dynVars_,iter);
				}
				//model_.doMeasurements(iter); FIXME
			}
		}
		
		const ParametersType params_;
		ModelType& model_;
		DynVarsType& dynVars_;
		ConcurrencyType& concurrency_;
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
