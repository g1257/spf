
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
#include "ProgressIndicator.h"
#include "MonteCarlo.h"

namespace Spf {
	
	template<typename ParametersType,typename AlgorithmType,typename ModelType,typename ConcurrencyType,typename RandomNumberGeneratorType>
	class Engine {
		
		typedef typename ParametersType::FieldType FieldType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef MonteCarlo<ParametersType,ModelType,AlgorithmType,RandomNumberGeneratorType> MonteCarloType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,AlgorithmType& algorithm,ConcurrencyType& concurrency) 
			: params_(params),model_(model),dynVars_(model.dynVars()),concurrency_(concurrency),fout_(params_.filename.c_str()),
				  progress_("Engine",concurrency.rank()),monteCarlo_(params,model,algorithm)
		{
		}
				
		void main()
		{
			thermalize();
			// announce thermalization done
			measure();
			// announce measurements done
			finalize();
		}
		
		private:
		
		void thermalize()
		{
			size_t acc = 0;
			for (size_t iter=0;iter<params_.iterTherm;iter++) {
				utils::printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				acc += doMonteCarlo(dynVars_,iter);
			}
			if (params_.iterTherm ==0) return;
			size_t pp = 100*acc/params_.iterTherm;
			std::string s = "Thermalization finished. " + utils::ttos(acc) + " or " + utils::ttos(pp) + "%";
			progress_.printline(s,fout_);
		}
		
		void measure()
		{
			size_t acc = 0;
			size_t counter = 0;
			for (size_t iter=0;iter<params_.iterEffective;iter++) {
				utils::printProgress(iter,params_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					acc += doMonteCarlo(dynVars_,iter);
					counter++;
				}
				model_.doMeasurements(dynVars_,iter,fout_);
				if (counter==0) continue;
				size_t pp = 100*acc/(counter*model_.totalFlips());
				std::string s = "Acceptance: " + utils::ttos(pp) + "%";
				progress_.printline(s,fout_);
			}
		}
		
		void finalize()
		{
			fout_<<model_;
			fout_<<dynVars_;
			std::vector<FieldType> eigs;
			model_.fillAndDiag(eigs,dynVars_);
			fout_<<"Eigenvalues\n";
			fout_<<eigs;
			std::cerr<<"\n";
		}
		
		size_t doMonteCarlo(DynVarsType& dynVars, size_t iter)
		{
			return monteCarlo_(dynVars,iter);
			
		}
		
		const ParametersType params_;
		ModelType& model_;
		DynVarsType& dynVars_;
		ConcurrencyType& concurrency_;
		std::ofstream fout_;
		ProgressIndicatorType progress_;
		MonteCarloType monteCarlo_;
		
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
