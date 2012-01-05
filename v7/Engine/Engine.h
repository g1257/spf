
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

namespace Spf {
	
	template<typename ParametersType,typename ModelType,typename IoInType,typename RngType>
	class Engine {
		
		typedef typename ParametersType::RealType RealType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::ConcurrencyType ConcurrencyType;
		typedef typename ConcurrencyType::CommType CommType;
		typedef Packer<RealType,PsimagLite::IoSimple::Out,ConcurrencyType> PackerType;
		typedef SaveConfigs<ParametersType,DynVarsType> SaveConfigsType;
		typedef Spf::GreenFunctionTpem<ParametersType,ModelType,RngType> GreenFunctionTpemType;
		typedef Spf::GreenFunctionDiag<ParametersType,ModelType,RngType> GreenFunctionDiagType;
		typedef AlgorithmFactory<GreenFunctionDiagType,GreenFunctionTpemType> AlgorithmFactoryType;

		static const std::string license_;

	public:

		Engine(ParametersType& params,
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
			const std::string opts = params.options;
			bool tpem = (opts.find("tpem")!=std::string::npos);
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

		static const std::string& license()
		{
			return license_;
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
				std::string algorithmicError = "DISABLED";
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
			std::string s(ctime(&t));
			ioOut_<<s;
			ioOut_<<params_;
			ioOut_<<model_;
		}

		void finalize()
		{
			model_.finalize(ioOut_,comm_.second);
			ioOut_<<"#FinalClassicalFieldConfiguration:\n";
			ioOut_<<dynVars_;
			bool printgf = (params_.options.find("printgf")!=std::string::npos);
			if (printgf) printGf();
			time_t t = time(0);
			std::string s(ctime(&t));
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

		void doMonteCarlo(std::vector<PairType>& accepted,DynVarsType& dynVars, size_t iter)
		{
			typedef typename DynVarsType::OperationsType0 OperationsType0;
			typedef typename DynVarsType::Type0 Type0;
			typedef MonteCarlo<ParametersType,OperationsType0,AlgorithmFactoryType,
			                   RngType,Type0> MonteCarloType0;
			AlgorithmFactoryType algorithm(gfDiag_,gfTpem_);

			MonteCarloType0 monteCarlo0(params_,model_.ops((OperationsType0*)0),algorithm,rng_);
			Type0& spinPart = dynVars.getField((Type0*)0);
			PairType res= monteCarlo0(spinPart,iter); //concurrency_,comm_.second);
			accepted[0].first += res.first;
			accepted[0].second += res.second;
			
			if (dynVars.size()==1) return;
			
			typedef typename DynVarsType::OperationsType1 OperationsType1;
			typedef typename DynVarsType::Type1 Type1;
			typedef MonteCarlo<ParametersType,OperationsType1,AlgorithmFactoryType,RngType,
   				Type1> MonteCarloType1;
			
			MonteCarloType1 monteCarlo1(params_,model_.ops((OperationsType1*)0),algorithm,rng_);
			Type1& phononPart = dynVars.getField((Type1*)0);
			res= monteCarlo1(phononPart,iter); //concurrency_,comm_.second);
			accepted[1].first += res.first;
			accepted[1].second += res.second;
			
		}

		void printProgress(const std::vector<PairType>& accepted,
		                   const std::string& algorithmicError = "DISABLED",
		                   PackerType* packer = 0)
		{
			for (size_t i=0;i<dynVars_.size();i++) {
				if (accepted[i].second==0) continue;
				std::string s1=  "Acceptance " + dynVars_.name(i) + "=";
				size_t pp = 100*accepted[i].first/accepted[i].second;
				std::string s2=  "AcceptancePercentage " + dynVars_.name(i) + "=%";
				std::string s3 = "AlgorithmicError=";

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
	
template<typename ParametersType,typename ModelType,typename IoInType,typename RngType>
const std::string Engine<ParametersType,ModelType,IoInType,RngType>::license_=
"Copyright (c) 2009-2011, UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[SPF, Version 7.0.0]\n"
"\n"
"*********************************************************\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. \n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"*********************************************************\n"
"\n"
"\n";
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
