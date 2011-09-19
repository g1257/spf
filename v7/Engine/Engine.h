
/** \ingroup SPF */
/*@{*/

/*! \file Engine.h
 *
 *  SPF Engine
 *
 */
#ifndef SPF_ENGINE_H
#define SPF_ENGINE_H
#include <fstream>
#include <iostream>
#include "IoSimple.h"
#include "ProgressIndicator.h" //in PsimagLite
#include "TypeToString.h" // in PsimagLite
#include "MonteCarlo.h"
#include "Packer.h"
#include "SaveConfigs.h"

namespace Spf {
	
	template<typename ParametersType,typename AlgorithmType,typename ModelType,
 		typename ConcurrencyType,typename RandomNumberGeneratorType,typename GreenFunctionType>
	class Engine {
		
		typedef typename ParametersType::FieldType FieldType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
		typedef std::pair<size_t,size_t> PairType;
		typedef Packer<FieldType,PsimagLite::IoSimple::Out,ConcurrencyType> PackerType;
		typedef SaveConfigs<ParametersType,DynVarsType> SaveConfigsType;
	public:
			
		Engine(ParametersType& params,ModelType& model,
		       AlgorithmType& algorithm,
		       ConcurrencyType& concurrency) 
		: params_(params),algorithm_(algorithm),model_(model),
		  dynVars_(model.dynVars()),concurrency_(concurrency),
		  ioOut_(params_.filename,concurrency_.rank()),progress_("Engine",concurrency.rank()),
		  rng_(params.randomSeed,concurrency_.rank(),concurrency_.nprocs()),saveConfigs_(params_,dynVars_,concurrency.rank())
		{
			size_t nprocs = concurrency_.nprocs();
			size_t temp = params_.iterEffective/nprocs;
			if (temp * nprocs != params_.iterEffective) {
				std::string s = "numberOfProcessors must be a divisor of params.iterEffective\n";
				std::runtime_error(s.c_str());
			}
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
				printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				doMonteCarlo(accepted,dynVars_,iter);
			}
			std::cerr<<"\n";
			printProgress(accepted);
		}
		
		void measure()
		{
			std::vector<std::pair<size_t,size_t> > accepted(dynVars_.size());
			concurrency_.loopCreate(params_.iterEffective);
			size_t iter=0;
			while(concurrency_.loop(iter)) {
				printProgress(iter,params_.iterEffective,10,'*',concurrency_.rank());
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					doMonteCarlo(accepted,dynVars_,iter);
				}
				GreenFunctionType greenFunction(params_,algorithm_,model_.hilbertSize());
				PackerType packer(ioOut_,concurrency_);
				model_.doMeasurements(greenFunction,iter,packer);
				saveConfigs_(iter); 
				printProgress(accepted,&packer);
			}
			std::cerr<<"\n";
		}
		
		void writeHeader()
		{
			ioOut_<<"#This is SPF v7\n";
			time_t t = time(0);
			std::string s(ctime(&t));
			ioOut_<<s;
			ioOut_<<params_;
			ioOut_<<model_;
		
		}
		
		void finalize()
		{
			model_.finalize(ioOut_);
			ioOut_<<"#FinalClassicalFieldConfiguration:\n";
			ioOut_<<dynVars_;
			ioOut_<<"#AlgorithmRelated:\n";
			ioOut_<<algorithm_;
			time_t t = time(0);
			std::string s(ctime(&t));
			ioOut_<<s;
			ioOut_<<"#EOF\n";
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
		
		void printProgress(const std::vector<PairType>& accepted,PackerType* packer = 0)
		{
			for (size_t i=0;i<dynVars_.size();i++) {
				if (accepted[i].second==0) continue;
				std::string s1=  "Acceptance " + dynVars_.name(i) + "=";
				size_t pp = 100*accepted[i].first/accepted[i].second;
				std::string s2=  "AcceptancePercentage " + dynVars_.name(i) + "=%";
				
				if (packer) {
					packer->pack(s1,accepted[i].first);
					packer->pack(s2,pp);
				} else {
					progress_.printline(s1+ttos(accepted[i].first),ioOut_);
					progress_.printline(s2+ttos(pp),ioOut_);
				}
			}
		}
		
		void printProgress(int i,int total,int nMarks,char mark,int option)
		{
			int every=total/nMarks;
			if (every<=0 || i<=0) return;
			if (i%every ==0) {
				std::cerr<<mark;
				std::cerr.flush();
			}
		}

		const ParametersType& params_;
		AlgorithmType& algorithm_;
		ModelType& model_;
		DynVarsType& dynVars_;
		ConcurrencyType& concurrency_;
		PsimagLite::IoSimple::Out ioOut_;
		ProgressIndicatorType progress_;
		RandomNumberGeneratorType rng_;
		SaveConfigsType saveConfigs_;
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
