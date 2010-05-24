
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
		typedef std::pair<size_t,size_t> PairType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,AlgorithmType& algorithm,ConcurrencyType& concurrency) 
			: params_(params),algorithm_(algorithm),model_(model),dynVars_(model.dynVars()),
				  concurrency_(concurrency),fout_(params_.filename.c_str()),
				  progress_("Engine",concurrency.rank())
		{
			rng_.seed(params_.randomSeed);
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
				GreenFunctionType greenFunction(params_,algorithm_,model_.hilbertSize());
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
			fout_<<params_;
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
			typedef typename DynVarsType::OperationsType0 OperationsType0;
			typedef typename DynVarsType::Type0 Type0;
			typedef MonteCarlo<ParametersType,OperationsType0,AlgorithmType,RandomNumberGeneratorType,
					Type0> MonteCarloType0;

			MonteCarloType0 monteCarlo0(params_,model_.ops((OperationsType0*)0),algorithm_,rng_);
			Type0& spinPart = dynVars.getField((Type0*)0);
			PairType res= monteCarlo0(spinPart,iter); // (accepted, totalflips)
			accepted[0].first += res.first;
			accepted[0].second += res.second;
			
			if (dynVars.size()==1) return;
			
			typedef typename DynVarsType::OperationsType1 OperationsType1;
			typedef typename DynVarsType::Type1 Type1;
			typedef MonteCarlo<ParametersType,OperationsType1,AlgorithmType,RandomNumberGeneratorType,
   				Type1> MonteCarloType1;
			
			MonteCarloType1 monteCarlo1(params_,model_.ops((OperationsType1*)0),algorithm_,rng_);
			Type1& phononPart = dynVars.getField((Type1*)0);
			res= monteCarlo1(phononPart,iter); // (accepted, totalflips)
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
		
		const ParametersType& params_;
		AlgorithmType& algorithm_;
		ModelType& model_;
		DynVarsType& dynVars_;
		ConcurrencyType& concurrency_;
		std::ofstream fout_;
		ProgressIndicatorType progress_;
		RandomNumberGeneratorType rng_;
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
