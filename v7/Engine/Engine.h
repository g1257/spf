
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
		typedef std::pair<size_t,size_t> PairType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,AlgorithmType& algorithm,ConcurrencyType& concurrency) 
			: params_(params),algorithm_(algorithm),model_(model),dynVars_(model.dynVars()),
				  concurrency_(concurrency),fout_(params_.filename.c_str()),
				  progress_("Engine",concurrency.rank()),monteCarlo_(params,model,algorithm)
		{
			writeHeader();
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
			std::vector<PairType> accepted(dynVars_.size());
			for (size_t iter=0;iter<params_.iterTherm;iter++) {
				utils::printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				doMonteCarlo(accepted,dynVars_,iter);
			}
			std::cerr<<"\n";
			printProgress(accepted);
			
// 			if (params_.iterTherm ==0) return;
// 			std::string s = "Thermalization finished. ";
// 			progress_.printline(s,fout_);
// 			for (size_t i=0;i<dynVars_.size();i++) {
// 				size_t pp = 100*accepted[i]/params_.iterTherm;
// 				s=  "Acceptance: " + dynVars_.name(i) + " " + utils::ttos(accepted[i]) +
// 							" or " + utils::ttos(pp) + "%";
// 				progress_.printline(s,fout_);
// 			}
		}
		
		void measure()
		{
			std::vector<std::pair<size_t,size_t> > accepted(dynVars_.size());
			for (size_t iter=0;iter<params_.iterEffective;iter++) {
				utils::printProgress(iter,params_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					doMonteCarlo(accepted,dynVars_,iter);
				}
				GreenFunctionType greenFunction(algorithm_);
				model_.doMeasurements(greenFunction,iter,fout_);
				printProgress(accepted);
			}
			std::cerr<<"\n";
		}
		
		void writeHeader()
		{
			fout_<<"#This is SPF v7\n";
			time_t t = time(0);
			fout_<<ctime(&t);
			fout_<<model_;
		}
		
		void finalize()
		{
			model_.finalize(fout_);
			fout_<<"#FinalClassicalFieldConfiguration:\n";
			fout_<<dynVars_;
			fout_<<"#AlgorithmRelated:\n";
			fout_<<algorithm_;
			time_t t = time(0);
			fout_<<ctime(&t);
			fout_<<"#EOF\n";
		}
		
		
		void doMonteCarlo(std::vector<PairType>& accepted,DynVarsType& dynVars, size_t iter)
		{
			
			PairType res= monteCarlo_(dynVars.template getField<0,typename DynVarsType::Type0>(),iter); // (accepted, totalflips)
			accepted[0].first += res.first;
			accepted[0].second += res.second;
			
			if (dynVars.size()==1) return;
			
			res= monteCarlo_(dynVars.template getField<1,typename DynVarsType::Type1>(),iter); // (accepted, totalflips)
			accepted[1].first += res.first;
			accepted[1].second += res.second;
			
		}
		
		void printProgress(const std::vector<PairType>& accepted)
		{
			for (size_t i=0;i<dynVars_.size();i++) {
				if (accepted[i].second==0) continue;
				size_t pp = 100*accepted[i].first/accepted[i].second;
				std::string s=  "Acceptance: " + dynVars_.name(i) + " " + utils::ttos(accepted[i].first) +
						" or " + utils::ttos(pp) + "%";
				progress_.printline(s,fout_);
			}
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
