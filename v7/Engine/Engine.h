
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
#include "Range.h" // in PsimagLite
#include "MonteCarlo.h"
#include "Packer.h"
#include "SaveConfigs.h"
#include "AlgorithmFactory.h"
#include "GreenFunctionTpem.h"
#include "GreenFunctionDiag.h"
#include <time.h>
#include "loki/Typelist.h"

namespace Spf {

	template<typename ParametersType,typename ModelType,typename IoInType,typename RngType>
	class Engine {
		
		typedef typename ParametersType::RealType RealType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef typename DynVarsType::OperationsList OperationsListType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::ConcurrencyType ConcurrencyType;
		typedef typename ConcurrencyType::CommType CommType;
		typedef Packer<RealType,PsimagLite::IoSimple::Out,ConcurrencyType> PackerType;
		typedef SaveConfigs<ParametersType,DynVarsType> SaveConfigsType;
		typedef Spf::GreenFunctionTpem<ParametersType,ModelType,RngType> GreenFunctionTpemType;
		typedef Spf::GreenFunctionDiag<ParametersType,ModelType,RngType> GreenFunctionDiagType;
		typedef AlgorithmFactory<GreenFunctionDiagType,GreenFunctionTpemType> AlgorithmFactoryType;

	public:

		Engine(const ParametersType& params,
		       ModelType& model,
		       IoInType& io,
		       ConcurrencyType& concurrency) 
		: params_(params),
		  model_(model),
		  concurrency_(concurrency),
		  comm_(concurrency_.newCommFromSegments(params_.coresForKernel)),
		  gfTpem_(0),
		  gfDiag_(0),
		  dynVars_(model_.dynVars()),
		  ioOut_(params_.filename,
		  concurrency_.rank()),
		  progress_("Engine",concurrency.rank()),
		  rng_(params.randomSeed,concurrency_.rank(comm_.second),concurrency_.nprocs(comm_.second)),
		  saveConfigs_(params_,dynVars_,concurrency.rank(comm_.first))
		{
			const PsimagLite::String opts = params.options;
			bool tpem = (opts.find("tpem")!=PsimagLite::String::npos);
			if (tpem) {
				gfTpem_ = new GreenFunctionTpemType(params,model,io,comm_.first);
			} else {
				gfDiag_ = new GreenFunctionDiagType(params,model,io,comm_.first);
			}
			writeHeader();
		}

		~Engine()
		{
			if (gfTpem_) delete gfTpem_;
			if (gfDiag_) delete gfDiag_;
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
			size_t fieldsToIntegrate = Loki::TL::Length<OperationsListType>::value;
			PsimagLite::Vector<PairType>::Type accepted(fieldsToIntegrate);
			for (size_t iter=0;iter<params_.iterTherm;iter++) {
				printProgress(iter,params_.iterTherm,10,'*',concurrency_.rank());
				doMonteCarlo(accepted,dynVars_,iter);
			}
			std::cerr<<"\n";
			printProgress(accepted);
		}

		void measure()
		{
			size_t fieldsToIntegrate = Loki::TL::Length<OperationsListType>::value;
			typedef std::pair<size_t,size_t> PairType;
			typename PsimagLite::Vector<PairType>::Type accepted(fieldsToIntegrate);

			bool isStrict = true;
			PsimagLite::Range<ConcurrencyType> range(0,params_.iterEffective,
			             concurrency_,comm_.second,isStrict);
						 
			for (;!range.end();range.next()) {
				size_t iter=range.index();
				printProgress(iter,params_.iterEffective,10,'*',
				         concurrency_.rank(comm_.second));
				for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
					doMonteCarlo(accepted,dynVars_,iter);
				}
				PackerType packer(ioOut_,concurrency_,comm_.second);
				PsimagLite::String algorithmicError = "DISABLED";
				if (gfDiag_) {
					gfDiag_->measure();
					model_.doMeasurements(*gfDiag_,iter,packer);
				} else {
					gfTpem_->measure();
					model_.doMeasurements(*gfTpem_,iter,packer);
					algorithmicError = gfTpem_->error();
				}
				saveConfigs_(iter); 
				printProgress(accepted,algorithmicError,&packer);
			}
			std::cerr<<"\n";
		}

		void writeHeader()
		{
			ioOut_<<"#This is SPF v7\n";
			time_t t = time(0);
			PsimagLite::String s(ctime(&t));
			ioOut_<<s;
			ioOut_<<params_;
			ioOut_<<model_;
		}

		void finalize()
		{
			model_.finalize(ioOut_,comm_.second);
			ioOut_<<"#FinalClassicalFieldConfiguration:\n";
			ioOut_<<dynVars_;
			bool printgf = (params_.options.find("printgf")!=PsimagLite::String::npos);
			if (printgf) printGf();
			time_t t = time(0);
			PsimagLite::String s(ctime(&t));
			ioOut_<<s;
			ioOut_<<"#EOF\n";
		}

		void printGf()
		{
			ioOut_<<"#AlgorithmRelated:\n";
			if (gfDiag_) {
				ioOut_<<(*gfDiag_);
			} else {
				ioOut_<<(*gfTpem_);
			}
		}

		void doMonteCarlo(PsimagLite::Vector<PairType>::Type& accepted,
		                  DynVarsType& dynVars, size_t iter)
		{
			typedef typename DynVarsType::OperationsList OperationsListType;
			AlgorithmFactoryType algorithm(gfDiag_,gfTpem_);

			MonteCarloLoop<RngType,
			        ParametersType,
			        ModelType,
			        AlgorithmFactoryType,
			        OperationsListType,
			        Loki::TL::Length<OperationsListType>::value-1>
			        ::loop(rng_,params_,algorithm,model_,dynVars,accepted,iter);
		}

		void printProgress(const PsimagLite::Vector<PairType>::Type& accepted,
		                   const PsimagLite::String& algorithmicError = "DISABLED",
		                   PackerType* packer = 0)
		{
			size_t fieldsToIntegrate = Loki::TL::Length<OperationsListType>::value;
			for (size_t i=0;i<fieldsToIntegrate;i++) {
				if (accepted[i].second==0) continue;
				PsimagLite::String s1=  "Acceptance " + dynVars_.name(i) + "=";
				size_t pp = 100*accepted[i].first/accepted[i].second;
				PsimagLite::String s2=  "AcceptancePercentage " + dynVars_.name(i) + "=%";
				PsimagLite::String s3 = "AlgorithmicError=";

				if (packer) {
					packer->pack(s1,accepted[i].first);
					packer->pack(s2,pp);
					if (algorithmicError!="DISABLED")
						packer->pack(s3,algorithmicError);
				} else {
					progress_.printline(s1+ttos(accepted[i].first),ioOut_);
					progress_.printline(s2+ttos(pp),ioOut_);
					if (algorithmicError!="DISABLED")
						progress_.printline(s3+algorithmicError,ioOut_);
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
		ModelType& model_;
		ConcurrencyType& concurrency_;
		/** comm_ is a pair of communicators that the Engine owns
		  * comm_.first is the communicator for the kernel, which
		  *   means either parallel TPEM or parallel diagonalization
		  * comm_.second is the communicator for the Monte Carlo
		  *  meaning that MC measurements are parallel over comm_.second
		  */ 
		std::pair<CommType,CommType> comm_;
		GreenFunctionTpemType* gfTpem_; // we own it, we new it, and we delete it
		GreenFunctionDiagType* gfDiag_; // we own it, we new it, and we delete it
		DynVarsType& dynVars_;
		PsimagLite::IoSimple::Out ioOut_;
		ProgressIndicatorType progress_;
		RngType rng_;
		SaveConfigsType saveConfigs_;
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
