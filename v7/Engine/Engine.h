
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
	
	template<typename ParametersType,typename AlgorithmType,typename ModelType,
 		typename ConcurrencyType,typename RandomNumberGeneratorType,typename GreenFunctionType>
	class Engine {
		
		typedef typename ParametersType::FieldType FieldType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef MonteCarlo<ParametersType,ModelType,AlgorithmType,RandomNumberGeneratorType> MonteCarloType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,AlgorithmType& algorithm,ConcurrencyType& concurrency) 
			: params_(params),algorithm_(algorithm),model_(model),dynVars_(model.dynVars()),
				  concurrency_(concurrency),fout_(params_.filename.c_str()),
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
			std::vector<size_t> accepted(dynVars_.size());
			for (size_t iter=0;iter<params_.iterTherm;iter++) {
				utils::printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				doMonteCarlo(accepted,dynVars_,iter);
			}
			if (params_.iterTherm ==0) return;
			std::string s = "Thermalization finished. ";
			progress_.printline(s,fout_);
			for (size_t i=0;i<dynVars_.size();i++) {
				size_t pp = 100*accepted[i]/params_.iterTherm;
				s=  "Acceptance: " + dynVars_.name(i) + " " + utils::ttos(accepted[i]) +
							" or " + utils::ttos(pp) + "%";
				progress_.printline(s,fout_);
			}
		}
		
		void measure()
		{
			std::vector<size_t> accepted(dynVars_.size());
			size_t counter = 0;
			for (size_t iter=0;iter<params_.iterEffective;iter++) {
				utils::printProgress(iter,params_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					doMonteCarlo(accepted,dynVars_,iter);
					counter++;
				}
				GreenFunctionType greenFunction(algorithm_);
				model_.doMeasurements(greenFunction,iter,fout_);
				if (counter==0) continue;
				for (size_t i=0;i<dynVars_.size();i++) {
					size_t pp = 100*accepted[i]/counter;
					std::string s=  "Acceptance: " + dynVars_.name(i) + " " + utils::ttos(accepted[i]) +
							" or " + utils::ttos(pp) + "%";
					progress_.printline(s,fout_);
				}
			}
		}
		
		void finalize()
		{
			fout_<<model_;
			fout_<<dynVars_;
			fout_<<algorithm_;
			std::cerr<<"\n";
		}
		
		void doMonteCarlo(std::vector<size_t> accepted,DynVarsType& dynVars, size_t iter)
		{
			for (size_t i=0;i<dynVars.size();i++) {
				accepted[i] += monteCarlo_(dynVars.getField(),iter);
			
		}
		
		const ParametersType params_;
		AlgorithmType& algorithm_;
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
