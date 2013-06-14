
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
#include "AlgorithmFactory.h"
#include "GreenFunctionTpem.h"
#include "GreenFunctionDiag.h"
#include <time.h>
#include "loki/Typelist.h"
#include "Concurrency.h"

namespace Spf {

	template<typename ParametersType,typename ModelType,typename IoInType,typename RngType>
	class Engine {
		
		typedef typename ParametersType::RealType RealType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef typename DynVarsType::OperationsList OperationsListType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
		typedef std::pair<size_t,size_t> PairType;
		typedef Packer<RealType,PsimagLite::IoSimple::Out> PackerType;
		typedef SaveConfigs<ParametersType,DynVarsType> SaveConfigsType;
		typedef Spf::GreenFunctionTpem<ParametersType,ModelType,RngType> GreenFunctionTpemType;
		typedef Spf::GreenFunctionDiag<ParametersType,ModelType,RngType> GreenFunctionDiagType;
		typedef AlgorithmFactory<GreenFunctionDiagType,GreenFunctionTpemType> AlgorithmFactoryType;

		static const SizeType CONCURRENCY_RANK = 0;

		class MyLoop {

			typedef PsimagLite::Concurrency ConcurrencyType;
			typedef std::pair<size_t,size_t> PairType;

		public:

			MyLoop(const ParametersType& params,
			       ModelType& model,
			       IoInType& io)
			: params_(params),
			  model_(model),
			  gfTpem_(0),
			  gfDiag_(0),
			  dynVars_(model_.dynVars()),
			  ioOut_(params_.filename,
			  CONCURRENCY_RANK),
			  progress_("Engine",PsimagLite::Concurrency::rank()),
			  rng_(params.randomSeed,PsimagLite::Concurrency::rank(comm_.second),PsimagLite::Concurrency::nprocs(comm_.second)),
			  saveConfigs_(params_,dynVars_,PsimagLite::Concurrency::rank(comm_.first)),
			  fieldsToIntegrate(Loki::TL::Length<OperationsListType>::value),
		      accepted(fieldsToIntegrate)
			{
				const PsimagLite::String opts = params.options;
				bool tpem = (opts.find("tpem")!=PsimagLite::String::npos);
				if (tpem) {
					gfTpem_ = new GreenFunctionTpemType(params,model,io);
				} else {
					gfDiag_ = new GreenFunctionDiagType(params,model,io);
				}
				writeHeader();
			}

			~MyLoop()
			{
				if (gfTpem_) delete gfTpem_;
				if (gfDiag_) delete gfDiag_;
			}

			void thread_function_(SizeType threadNum,
			                      SizeType blockSize,
			                      SizeType total,
			                      typename ConcurrencyType::MutexType* myMutex)
			{
				for (SizeType p=0;p<blockSize;p++) {
					SizeType taskNumber = threadNum*blockSize + p;
					if (taskNumber>=total) break;

					size_t iter = taskNumber;
					printProgress(iter,params_.iterEffective,10,'*',
					         PsimagLite::Concurrency::rank(comm_.second));
					for (size_t iter2=0;iter2<params_.iterUnmeasured;iter2++) {
						doMonteCarlo(accepted,iter);
					}
					PackerType packer(ioOut_);
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
			}

			template<typename SomeParallelType>
			void gather(SomeParallelType& p)
			{
				if (ConcurrencyType::mode == ConcurrencyType::MPI) {
					p.allGather(accepted);
				}
			}

			void doMonteCarlo(PsimagLite::Vector<PairType>::Type& accepted,
			                  size_t iter)
			{
				typedef typename DynVarsType::OperationsList OperationsListType;
				AlgorithmFactoryType algorithm(gfDiag_,gfTpem_);

				MonteCarloLoop<RngType,
				        ParametersType,
				        ModelType,
				        AlgorithmFactoryType,
				        OperationsListType,
				        Loki::TL::Length<OperationsListType>::value-1>
				        ::loop(rng_,params_,algorithm,model_,dynVars_,accepted,iter);
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
				model_.finalize(ioOut_);
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

		private:

			const ParametersType& params_;
			ModelType& model_;
			/** comm_ is a pair of communicators that the Engine owns
			  * comm_.first is the communicator for the kernel, which
			  *   means either parallel TPEM or parallel diagonalization
			  * comm_.second is the communicator for the Monte Carlo
			  *  meaning that MC measurements are parallel over comm_.second
			  */
			std::pair<int,int> comm_;
			GreenFunctionTpemType* gfTpem_; // we own it, we new it, and we delete it
			GreenFunctionDiagType* gfDiag_; // we own it, we new it, and we delete it
			DynVarsType& dynVars_;
			PsimagLite::IoSimple::Out ioOut_;
			ProgressIndicatorType progress_;
			RngType rng_;
			SaveConfigsType saveConfigs_;
			size_t fieldsToIntegrate;
			typename PsimagLite::Vector<PairType>::Type accepted;
		};

	public:

		Engine(const ParametersType& params,
		       ModelType& model,
		       IoInType& io)
		: params_(params),
		  model_(model),
		  helper_(params,model,io)
		{}

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
				helper_.printProgress(iter,params_.iterTherm,10,'*',CONCURRENCY_RANK);
				helper_.doMonteCarlo(accepted,iter);
			}
			std::cerr<<"\n";
			helper_.printProgress(accepted);
		}

		void measure()
		{
			typedef PsimagLite::Parallelizer<MyLoop> ParallelizerType;
			ParallelizerType threadObject;

			ParallelizerType::setThreads(params_.npthreads);

			threadObject.loopCreate(params_.iterEffective,helper_);

			helper_.gather(threadObject);

			std::cerr<<"\n";
		}

		void writeHeader()
		{
			helper_.writeHeader();
		}

		void finalize()
		{
			helper_.finalize();
		}

		void printGf()
		{
			helper_.printGf();
		}

		const ParametersType& params_;
		ModelType& model_;
		MyLoop helper_;
	}; // Engine
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H
