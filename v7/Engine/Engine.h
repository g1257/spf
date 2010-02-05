
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
	
	template<typename ModelType,typename ConcurrencyType>
	class Engine {
	
		public:
			
		Engine(ModelType& model,ConcurrencyType& concurrency) : model_(model),concurrency_(concurrency)
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
			for (size_t iter=0;iter<mp_.iterTherm;iter++) {
				utils::printProgress(iter,mp_.iterTherm,10,'*',concurrency_.rank());
				model_.doMonteCarlo(iter);
			}
		}
		
		void measure()
		{
			for (iter=0;iter<ether.iterEffective;iter++) {
				utils::printProgress(iter,mp_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<mp_.iterUnmeasured;iter2++) {
					model_.doMonteCarlo(iter);
				}
				//model_.doMeasurements(iter); FIXME
			}
		}
		
		
	}; // Engine
} // namespace Spf

#endif // SPF_ENGINE_H
