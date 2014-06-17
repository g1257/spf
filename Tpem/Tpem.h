
/** \ingroup TPEM */
/*@{*/

/*! \file Tpem.h
 *
 *  
 *
 */
#ifndef TPEM_H
#define TPEM_H

#include "CrsMatrix.h"
#include "TpemSubspace.h"
#include "TpemFunctors.h"
#include "TypeToString.h"
#include <cassert>
#include "GslWrapper.h"
#include "ChebyshevFunctionExplicit.h"
#include <cmath>
#include "BLAS.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Tpem {

	template<typename TpemParametersType_,typename RealOrComplexType>
	class Tpem {

	public:

		typedef TpemParametersType_ TpemParametersType;
		typedef typename TpemParametersType::RealType RealType;
		typedef BaseFunctor<TpemParametersType> BaseFunctorType;
		typedef PsimagLite::CrsMatrix<RealOrComplexType> TpemSparseType;
		typedef TpemSubspace<RealType,TpemSparseType> TpemSubspaceType;
		typedef ActionFunctor<TpemParametersType> ActionFunctorType;
		typedef EnergyFunctor<TpemParametersType> EnergyFunctorType;
		typedef NumberFunctor<TpemParametersType> NumberFunctorType;
		typedef PsimagLite::GslWrapper GslWrapperType;
		// choose how to compute ChebyshevFunction below
		typedef PsimagLite::ChebyshevFunctionExplicit<RealType> ChebyshevFunctionType;
		//typedef PsimagLite::ChebyshevFunctionCached<RealType> ChebyshevFunctionType;
		//typedef PsimagLite::ChebyshevFunction<RealType> ChebyshevFunctionType;
		
		enum {NO_VERBOSE,YES_VERBOSE};
		static const SizeType verbose_ = NO_VERBOSE;

		struct MyFunctionParams {
			MyFunctionParams(const BaseFunctorType& functor1)
			: functor(functor1) { }

			const BaseFunctorType& functor;
			SizeType m;
		};

		typedef MyFunctionParams MyFunctionParamsType;

		class MyLoop {

			typedef PsimagLite::Concurrency ConcurrencyType;

		public:

			MyLoop(std::vector<RealType> &vobs1,
			       const BaseFunctorType& obsFunc,
			       const GslWrapperType& gslWrapper)
			    : vobs(vobs1),
			      gslWrapper_(gslWrapper),
			      pts(2),
			      epsabs(1e-9),
			      epsrel(1e-9),
			      limit(1000000),
			      workspace(gslWrapper_.gsl_integration_workspace_alloc(limit+2)),
			      result(0),
			      abserr(0),
			      params(obsFunc)
			{
				pts[0]= -1.0;
				pts[1] = 1.0;

				f.function= &Tpem<TpemParametersType,RealOrComplexType>::myFunction;

				f.params = &params;
			}

			~MyLoop()
			{
				gslWrapper_.gsl_integration_workspace_free (workspace);
			}

			void thread_function_(SizeType threadNum,
			                      SizeType blockSize,
			                      SizeType total,
			                      typename ConcurrencyType::MutexType* myMutex)
			{
				SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
				SizeType npthreads = ConcurrencyType::npthreads;
				for (SizeType p=0;p<blockSize;p++) {
					SizeType taskNumber = (threadNum+npthreads*mpiRank)*blockSize + p;
					if (taskNumber>=total) break;

					//std::cout<<"This is thread number "<<threadNum;
					//std::cout<<" and taskNumber="<<taskNumber<<"\n";

					params.m = taskNumber;
					gslWrapper_.gsl_integration_qagp(&f,&(pts[0]),pts.size(),epsabs,epsrel,limit,workspace,&result,&abserr);
					//gsl_integration_qag(&f,pts[0],pts[1],epsabs,epsrel,limit,key,workspace,&result,&abserr);
					//gsl_integration_qags(&f,pts[0],pts[1],epsabs,epsrel,limit,workspace,&result,&abserr);

					if (std::isinf(result) || std::isnan(result)) {
						vobs[params.m] = 0;
						continue;
					}
					vobs[params.m] = result;
				}
			}

			void gather()
			{
				PsimagLite::MPI::pointByPointGather(vobs);
			}

		private:

			typename PsimagLite::Vector<RealType>::Type &vobs;
			const GslWrapperType& gslWrapper_;
			typename PsimagLite::Vector<RealType>::Type pts;
			RealType epsabs;
			RealType epsrel;
			SizeType limit;
			GslWrapperType::gsl_integration_workspace *workspace;
			RealType result;
			RealType abserr;
			GslWrapperType::gsl_function f;
			MyFunctionParamsType params;
		};

		class MyLoop2 {

			typedef PsimagLite::Concurrency ConcurrencyType;

		public:

			MyLoop2(TpemSparseType& matrix1,
			        std::vector<RealType>& moment1,
			        const TpemParametersType& tpemParameters)
			    : matrix(matrix1),moment(moment1),tpemParameters_(tpemParameters)
			{
				SizeType n=moment.size();

				for (SizeType i = 0; i < n; i++) moment[i] = 0.0;
				assert(matrix.row()==matrix.col());
				moment[0] = matrix.row();
			}

			void thread_function_(SizeType threadNum,
			                      SizeType blockSize,
			                      SizeType total,
			                      typename ConcurrencyType::MutexType* myMutex)
			{
				SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
				SizeType npthreads = ConcurrencyType::npthreads;
				for (SizeType p=0;p<blockSize;p++) {
					SizeType taskNumber = (threadNum+npthreads*mpiRank)*blockSize + p;
					if (taskNumber>=total) break;

					//std::cout<<"This is thread number "<<threadNum;
					//std::cout<<" and taskNumber="<<taskNumber<<"\n";

					diagonalElement(matrix, moment, taskNumber, tpemParameters_);
				}
			}

			void gather()
			{
				PsimagLite::MPI::pointByPointGather(moment);

				moment[0] = matrix.row();

				SizeType n=moment.size();

				for (SizeType i = 2; i < n; i += 2)
					moment[i] = 2.0 * moment[i] - moment[0];

				for (SizeType i = 3; i < n - 1; i += 2)
					moment[i] = 2.0 * moment[i] - moment[1];
			}

		private:

			TpemSparseType& matrix;
			typename PsimagLite::Vector<RealType>::Type& moment;
			const TpemParametersType& tpemParameters_;
		};

		class MyLoop3 {

			typedef PsimagLite::Concurrency ConcurrencyType;

		public:

			MyLoop3(const TpemSparseType& matrix0,
			        std::vector<RealType>& moment0,
			        const TpemSparseType& matrix1,
			        std::vector<RealType>& moment1,
			        TpemSubspaceType& info,
			        const TpemParametersType& tpemParameters)
			    : matrix0_(matrix0),
			      moment0_(moment0),
			      matrix1_(matrix1),
			      moment1_(moment1),
			      info_(info),
			      tpemParameters_(tpemParameters)
			{

			}

			void thread_function_(SizeType threadNum,
			                      SizeType blockSize,
			                      SizeType total,
			                      typename ConcurrencyType::MutexType* myMutex)
			{
				SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
				SizeType npthreads = ConcurrencyType::npthreads;
				for (SizeType p=0;p<blockSize;p++) {
					SizeType taskNumber = (threadNum+npthreads*mpiRank)*blockSize + p;
					if (taskNumber>=total) break;

					//std::cout<<"This is thread number "<<threadNum;
					//std::cout<<" and taskNumber="<<taskNumber<<"\n";

					SizeType ket = info_(taskNumber);
					diagonalElement(matrix0_, moment0_, ket,tpemParameters_);
					diagonalElement(matrix1_, moment1_, ket,tpemParameters_);
				}
			}

			void gather()
			{
				PsimagLite::MPI::pointByPointGather(moment0_);
				PsimagLite::MPI::pointByPointGather(moment1_);
			}

		private:

			const TpemSparseType& matrix0_;
			std::vector<RealType>& moment0_;
			const TpemSparseType& matrix1_;
			std::vector<RealType>& moment1_;
			TpemSubspaceType& info_;
			const TpemParametersType& tpemParameters_;
		};

		Tpem(const TpemParametersType& tpemParameters)
		    : tpemParameters_(tpemParameters)
		{
			gslWrapper_.gsl_set_error_handler(&my_handler);
		}

		void calcCoeffs(std::vector<RealType> &vobs,
		                const BaseFunctorType& obsFunc) const
		{
			typedef MyLoop HelperType;
			typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
			ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
			                              PsimagLite::MPI::COMM_WORLD);

			HelperType helper(vobs,obsFunc,gslWrapper_);

			threadObject.loopCreate(tpemParameters_.cutoff,helper);

			helper.gather();
		}

		void calcMoments(TpemSparseType& matrix,
		                 std::vector<RealType>& moment) const
		{
			typedef MyLoop2 Helper2Type;
			typedef PsimagLite::Parallelizer<Helper2Type> Parallelizer2Type;
			Parallelizer2Type threadObject(PsimagLite::Concurrency::npthreads,
			        PsimagLite::MPI::COMM_WORLD);

			Helper2Type helper(matrix,moment,tpemParameters_);

			threadObject.loopCreate(matrix.row(),helper);

			helper.gather();
		}

		void calcMomentsDiff(std::vector<RealType> &moments,
		                     const TpemSparseType& matrix0,
		                     const TpemSparseType& matrix1) const
		{
			assert(matrix0.row()==matrix0.col());
			TpemSubspaceType info(matrix0.row());
			SizeType n=moments.size();
			std::vector<RealType>  moment0(n,0.0), moment1(n,0.0);

			if (tpemParameters_.algorithm==TpemParametersType::TPEM) {
				subspaceForTrace(info,matrix0, matrix1, moment0, moment1,tpemParameters_);
				// reset all moments to zero because the subspace for trace changes them
				for (SizeType i = 0; i < n; i++) moment0[i] = moment1[i] = 0.0;
			} else {
				info.fill();
			}

			assert(matrix0.row()==matrix0.col());
			moment0[0] = moment1[0] =  matrix0.row();

			typedef MyLoop3 Helper3Type;
			typedef PsimagLite::Parallelizer<Helper3Type> Parallelizer3Type;
			Parallelizer3Type threadObject(PsimagLite::Concurrency::npthreads,
			                               PsimagLite::MPI::COMM_WORLD);

			Helper3Type helper(matrix0,moment0,matrix1,moment1,info,tpemParameters_);

			threadObject.loopCreate(info.top(),helper);

			helper.gather();

			moment0[0] = moment1[0] = matrix0.row();

			for (SizeType i = 2; i < n; i += 2) {
				moment0[i] = 2.0 * moment0[i] - moment0[0];
				moment1[i] = 2.0 * moment1[i] - moment1[0];
			}

			for (SizeType i = 3; i < n - 1; i += 2) {
				moment0[i] = 2.0 * moment0[i] - moment0[1];
				moment1[i] = 2.0 * moment1[i] - moment1[1];
			}

			for (SizeType i = 0; i < n; i++)
				moments[i] = moment0[i] - moment1[i];
		}

		RealType expand(const std::vector<RealType>& moments,
		                const std::vector<RealType>& coeffs,
		                SizeType start=0,
		                SizeType progressiveCutoff=0) const
		{
			assert(moments.size()==coeffs.size());

			if (progressiveCutoff==0) progressiveCutoff = coeffs.size();
			assert(progressiveCutoff<=coeffs.size());
			if (moments.size()<progressiveCutoff) progressiveCutoff = moments.size();
			if (start>=progressiveCutoff) start=0;
			double* dx = (double *) &(moments[start]);
			double* dy = (double *) &(coeffs[start]);
			RealType ret = psimag::BLAS::DOT(progressiveCutoff-start,dx,1,dy,1);
			/* for (SizeType i = 0; i < progressiveCutoff; ++i)
				ret += moments[i] * coeffs[i];
			*/
			assert(!std::isinf(ret) && !std::isnan(ret));
//			concurrency_.broadcast(ret,comm_);
			return ret;
		}

		const TpemParametersType& tpemParameters() const { return tpemParameters_; }

		static RealType myFunction (RealType x, void * p) 
		{
			MyFunctionParamsType* params = (MyFunctionParamsType*)p;

			/* return  params->funk(params->m,x);*/
			RealType tmp;
			RealType tmp2= (RealType)1.0/M_PI;
			int m=params->m;
			RealType factorAlpha = (m==0) ? 1 : 2;
			tmp = params->functor(x) *  tmp2 * factorAlpha * chebyshev_(m,x)/sqrt(1.0-x*x);
			return tmp;
		}

		static void my_handler (const char * reason, const char * file, int line, int gsl_errno)
		{
			if (verbose_==NO_VERBOSE) return;
			std::string s("GSL error handler called with reason=");
			s += std::string(reason) + std::string(" gsl_errno=")+ttos(gsl_errno);
			std::cerr<<s<<"\n";
		}

	private:

		static void diagonalElement(const TpemSparseType& matrix,
		                            std::vector<RealType> &moment,
		                            SizeType ket,
		                            const TpemParametersType& tpemParameters)
		{
			TpemSubspaceType work(matrix.row());
			if (tpemParameters.algorithm == TpemParametersType::TPEM) {
				diagonalElementTpem(matrix,moment,ket,work,tpemParameters.epsForProduct);
			} else if (tpemParameters.algorithm == TpemParametersType::PEM) {
				diagonalElementPem(matrix,moment,ket);
			} else {
				std::string s("tpem_diagonal_element: Unknown type: ");
				s += ttos(tpemParameters.algorithm) + "\n";
				throw std::runtime_error(s.c_str());
			}
		}

		static void diagonalElement(const TpemSparseType& matrix,
		                            std::vector<RealType> &moment,
		                            SizeType ket,
		                            TpemSubspaceType& info,
		                            const TpemParametersType& tpemParameters)
		{
			if (tpemParameters.algorithm != TpemParametersType::TPEM) {
				std::string s("tpem_diagonal_element: Expected algorithm==tpem");
				s += std::string("but found tpemType==");
				s += ttos(tpemParameters.algorithm) + "\n";
				throw std::runtime_error(s.c_str());
			}
			diagonalElementTpem(matrix,moment,ket,info,tpemParameters.epsForTrace);
		}

		static void diagonalElementTpem(const TpemSparseType& matrix,
		                                std::vector<RealType> &moment,
		                                SizeType ket,
		                                TpemSubspaceType& info,
		                                const RealType& eps)
		{
			assert(matrix.row()==matrix.col());
			std::vector<RealOrComplexType> tmp(matrix.row(),0.0);
			std::vector<RealOrComplexType> jm0(matrix.row(),0.0);
			std::vector<RealOrComplexType> jm1(matrix.row(),0.0);

			jm0[ket] = 1.0;/* set |j,0> */

			info.clear();
			info.push(ket);

			/* calculate |j,1> = X|j,0> */
			info.sparseProduct(matrix, jm1, jm0,eps);

			RealOrComplexType sum1 = std::conj(jm0[ket]) *  jm1[ket];
			moment[1] += std::real(sum1);

			RealOrComplexType sum2 = 0.0;
			for (SizeType p = 0; p < info.top(); p++)
				sum2 += std::conj(jm1[info(p)]) * jm1[info(p)];
			moment[2] += std::real (sum2);

			/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
			* 
			* begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
			* end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
			* ...
			* begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
			* end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
			*/
			SizeType n=moment.size();
			for (SizeType m = 2; m < n / 2; m++) {
				/* calculate |tmp> = X|jm1> */
				info.sparseProduct(matrix,tmp, jm1,eps);
				RealOrComplexType sum1 = 0.0;
				RealOrComplexType sum2 = 0.0;
				for (SizeType p = 0; p < info.top(); p++) {
					SizeType i = info(p);
					RealOrComplexType keep = tmp[i] + tmp[i]  - jm0[i];
					/* for moment[2 * m    ] */
					sum1 += std::conj(keep) * keep;
					/* for moment[2 * m - 1] */
					sum2 += std::conj(keep) * jm1[i];
					/* set |j,m-1> and |j,m-2> for next iteration */
					jm0[i] = jm1[i];
					jm1[i] = keep;
				}
				moment[m + m    ] += std::real (sum1);
				moment[m + m - 1] += std::real (sum2);
				assert(fabs(std::imag(sum2))<1e-6);
			}
		}

		static void  diagonalElementPem(const TpemSparseType& matrix,
		                                std::vector<RealType> &moment,
		                                SizeType ket)
		{
			assert(matrix.row()==matrix.col());
			std::vector<RealOrComplexType> tmp(matrix.row(),0.0);
			std::vector<RealOrComplexType> jm0(matrix.row(),0.0);
			std::vector<RealOrComplexType> jm1(matrix.row(),0.0);

			jm0[ket] = 1.0;	/* set |j,0> */

			/* calculate |j,1> = X|j,0> */
			TpemSubspaceType::sparseProductPem(matrix, jm1, jm0);

			RealOrComplexType sum1 = std::conj(jm0[ket]) * jm1[ket];
			moment[1] += std::real (sum1);

			RealOrComplexType sum2 =  0.0;
			for (SizeType j=0;j<matrix.row();j++)
				sum2 += std::conj(jm1[j]) * jm1[j];
			moment[2] += std::real (sum2);

			/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
			* 
			* begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
			* end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
			* ...
			* begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
			* end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
			*/
			SizeType n=moment.size();
			for (SizeType m = 2; m < n / 2; m++) {
				/* calculate |tmp> = X|jm1> */
				TpemSubspaceType::sparseProductPem(matrix, tmp, jm1);
				RealOrComplexType sum1 = 0.0;
				RealOrComplexType sum2 = 0.0;
				for (SizeType i=0;i<matrix.row();i++) {
					RealOrComplexType keep = tmp[i] + tmp[i] - jm0[i];
					/* for moment[2 * m    ] */
					sum1 +=  std::conj(keep) * keep;
					/* for moment[2 * m - 1] */
					sum2 += std::conj(keep) * jm1[i];
					/* set |j,m-1> and |j,m-2> for next iteration */
					jm0[i] = jm1[i];
					jm1[i] = keep;
				}
				moment[m + m    ] += std::real (sum1);
				moment[m + m - 1] += std::real (sum2);
			}
		}

		static void subspaceForTrace(TpemSubspaceType& info,
		                      const TpemSparseType& matrix0,
		                      const TpemSparseType& matrix1,
		                      std::vector<RealType>& moment0,
		                      std::vector<RealType>& moment1,
		                      const TpemParametersType& tpemParameters)
		{
			const std::vector<SizeType>& support = tpemParameters.support;
			info.clear();

			assert(matrix0.row()==matrix0.col());
			TpemSubspaceType work(matrix0.row());
			for (SizeType i = 0; i < support.size(); i++) {
				SizeType j = support[i];
				diagonalElement(matrix0, moment0, j,work,tpemParameters);
				for (SizeType p = 0; p < work.top(); p++) info.push(work(p));
				diagonalElement (matrix1, moment1, j,work,tpemParameters);
				for (SizeType p = 0; p < work.top(); p++) info.push(work(p));
			}
		}

		const TpemParametersType& tpemParameters_; // not the owner, just a ref
		GslWrapperType gslWrapper_;
		static ChebyshevFunctionType chebyshev_;
	}; // class Tpem

	template<typename TpemParametersType,typename RealOrComplexType>
	typename Tpem<TpemParametersType,RealOrComplexType>::ChebyshevFunctionType
	Tpem<TpemParametersType,RealOrComplexType>::chebyshev_;
} // namespace Tpem

/*@}*/
#endif // TPEM_H
