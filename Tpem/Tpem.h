
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
#include "Range.h"
#include <cmath>
#include "BLAS.h"

namespace Tpem {

	template<typename TpemParametersType_,typename RealOrComplexType,typename ConcurrencyType>
	class Tpem {

		typedef typename ConcurrencyType::CommType CommType;

		static const CommType COMM_WORLD;

	public:

		typedef PsimagLite::Range<ConcurrencyType> RangeType;
		typedef TpemParametersType_ TpemParametersType;
		typedef typename TpemParametersType::RealType RealType;
		typedef BaseFunctor<TpemParametersType> BaseFunctorType;
		typedef PsimagLite::CrsMatrix<RealOrComplexType> TpemSparseType;
		typedef TpemSubspace<RealType,TpemSparseType> TpemSubspaceType;
		typedef ActionFunctor<TpemParametersType> ActionFunctorType;
		typedef EnergyFunctor<TpemParametersType> EnergyFunctorType;
		typedef NumberFunctor<TpemParametersType> NumberFunctorType;
		typedef PsimagLite::GslWrapper GslWrapperType;
		
		enum {NO_VERBOSE,YES_VERBOSE};
		static const size_t verbose_ = NO_VERBOSE;
		
		struct MyFunctionParams {
			MyFunctionParams(const BaseFunctorType& functor1)
			: functor(functor1) { }

			const BaseFunctorType& functor;
			size_t m;
		};

		typedef MyFunctionParams MyFunctionParamsType;

		Tpem(const TpemParametersType& tpemParameters,
		     ConcurrencyType& concurrency,
			 CommType comm = COMM_WORLD)
		: tpemParameters_(tpemParameters),concurrency_(concurrency),comm_(comm)
		{
			gslWrapper_.gsl_set_error_handler(&my_handler);
		}

		void calcCoeffs(std::vector<RealType> &vobs,
		                const BaseFunctorType& obsFunc) const
		{
			std::vector<RealType> pts(2);
			pts[0]= -1.0;
			pts[1] = 1.0;
			//size_t npts = pts.size();
			RealType epsabs=1e-9;
			RealType epsrel=1e-9;

			size_t limit = 1000000;
			GslWrapperType::gsl_integration_workspace *workspace = 
			                   gslWrapper_.gsl_integration_workspace_alloc(limit+2);

			RealType result = 0,abserr = 0;

			GslWrapperType::gsl_function f;
			f.function= &Tpem<TpemParametersType,RealOrComplexType,ConcurrencyType>::myFunction;
			MyFunctionParamsType params(obsFunc);
			f.params = &params;
			//int key = GSL_INTEG_GAUSS61;
			RangeType range(0,tpemParameters_.cutoff,concurrency_,comm_);
			for (;!range.end();range.next()) {
				params.m = range.index();
				gslWrapper_.gsl_integration_qagp(&f,&(pts[0]),pts.size(),epsabs,epsrel,limit,workspace,&result,&abserr);
				//gsl_integration_qag(&f,pts[0],pts[1],epsabs,epsrel,limit,key,workspace,&result,&abserr);
				//gsl_integration_qags(&f,pts[0],pts[1],epsabs,epsrel,limit,workspace,&result,&abserr);
				
				if (std::isinf(result) || std::isnan(result)) {
					vobs[params.m] = 0;
					continue;
				}
				vobs[params.m] = result;
			}
			concurrency_.reduce(vobs,comm_);
			gslWrapper_.gsl_integration_workspace_free (workspace);
		}

		void calcMoments(TpemSparseType& matrix,
		                 std::vector<RealType>& moment) const
		{	
			size_t n=moment.size();
			std::vector<RealType> buf(n,0);

			for (size_t i = 0; i < n; i++) moment[i] = 0.0;
			moment[0] = matrix.rank();
			
			RangeType range(0,matrix.rank(),concurrency_,comm_);
			for (;!range.end();range.next()) {
				size_t i = range.index();
				diagonalElement(matrix, moment, i);
			}
			concurrency_.reduce(moment,comm_);
			moment[0] = matrix.rank();

			for (size_t i = 2; i < n; i += 2)
				moment[i] = 2.0 * moment[i] - moment[0];
				
			for (size_t i = 3; i < n - 1; i += 2)
				moment[i] = 2.0 * moment[i] - moment[1];
		}

		void calcMomentsDiff(std::vector<RealType> &moments,
		                     const TpemSparseType& matrix0,
		                     const TpemSparseType& matrix1) const
		{	
			TpemSubspaceType info(matrix0.rank());
			size_t n=moments.size();
			std::vector<RealType>  moment0(n,0.0), moment1(n,0.0);

			if (tpemParameters_.algorithm==TpemParametersType::TPEM) {
				subspaceForTrace(info,matrix0, matrix1, moment0, moment1);
				// reset all moments to zero because the subspace for trace changes them
				for (size_t i = 0; i < n; i++) moment0[i] = moment1[i] = 0.0;
			} else {
				info.fill();
			}

			moment0[0] = moment1[0] =  matrix0.rank();

			RangeType range(0,info.top(),concurrency_,comm_);
			// FIXME: we could use twice as many procs here
			for (;!range.end();range.next()) {
				size_t p= info(range.index());
				diagonalElement(matrix0, moment0, p);
				diagonalElement(matrix1, moment1, p);
			}
			concurrency_.reduce(moment0,comm_);
			concurrency_.reduce(moment1,comm_);
			moment0[0] = moment1[0] = matrix0.rank();

			for (size_t i = 2; i < n; i += 2) {
				moment0[i] = 2.0 * moment0[i] - moment0[0];
				moment1[i] = 2.0 * moment1[i] - moment1[0];
			}

			for (size_t i = 3; i < n - 1; i += 2) {
				moment0[i] = 2.0 * moment0[i] - moment0[1];
				moment1[i] = 2.0 * moment1[i] - moment1[1];
			}

			for (size_t i = 0; i < n; i++)
				moments[i] = moment0[i] - moment1[i];
		}

		RealType expand(const std::vector<RealType>& moments,
		                const std::vector<RealType>& coeffs,
						size_t progressiveCutoff=0) const
		{
			assert(moments.size()==coeffs.size());

			if (progressiveCutoff==0) progressiveCutoff = coeffs.size();
			assert(progressiveCutoff<=coeffs.size());
			if (moments.size()<progressiveCutoff) progressiveCutoff = moments.size();
			double* dx = (double *) &(moments[0]);
			double* dy = (double *) &(coeffs[0]);
			RealType ret = psimag::BLAS::DOT(progressiveCutoff,dx,1,dy,1);
			/* for (size_t i = 0; i < progressiveCutoff; ++i)
				ret += moments[i] * coeffs[i];
			*/
			assert(!std::isinf(ret) && !std::isnan(ret));
			concurrency_.broadcast(ret,comm_);
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
			tmp = params->functor(x) *  tmp2 * factorAlpha * chebyshev(m,x)/sqrt(1.0-x*x);
			return tmp;
		}

		static void my_handler (const char * reason, const char * file, int line, int gsl_errno)
		{
			if (verbose_==NO_VERBOSE) return;
			std::string s("GSL error handler called with reason=");
			s += std::string(reason) + std::string(" gsl_errno=")+ttos(gsl_errno);
			std::cerr<<s<<"\n";
		}

// 		ConcurrencyType& concurrency() { return concurrency_; }
// 
// 		CommType comm() const { return comm_; }

	private:

		static RealType chebyshev(int m,RealType x)
		{
			if (m==0) return 1;

			if (m==1) return x;
			
			if (m&1) {
				int p=(m-1)/2;
				return (2*chebyshev(p,x)*chebyshev(p+1,x)-x);
			}

			int pp = m/2;
			RealType tmp=chebyshev(pp,x);
			return (2*tmp*tmp-1);
		}

		void diagonalElement(const TpemSparseType& matrix,
		                     std::vector<RealType> &moment,
		                     size_t ket) const
		{
			TpemSubspaceType work(matrix.rank());
			if (tpemParameters_.algorithm == TpemParametersType::TPEM) {
				diagonalElementTpem(matrix,moment,ket,work,tpemParameters_.epsForProduct);
			} else if (tpemParameters_.algorithm == TpemParametersType::PEM) {
				diagonalElementPem(matrix,moment,ket);
			} else {
				std::string s("tpem_diagonal_element: Unknown type: ");
				s += ttos(tpemParameters_.algorithm) + "\n";
				throw std::runtime_error(s.c_str());
			}
		}

		void diagonalElement(const TpemSparseType& matrix,
		                     std::vector<RealType> &moment,
		                     size_t ket,
		                     TpemSubspaceType& info) const
		{
			if (tpemParameters_.algorithm != TpemParametersType::TPEM) {
				std::string s("tpem_diagonal_element: Expected algorithm==tpem");
				s += std::string("but found tpemType==");
				s += ttos(tpemParameters_.algorithm) + "\n";
				throw std::runtime_error(s.c_str());
			}
			diagonalElementTpem(matrix,moment,ket,info,tpemParameters_.epsForTrace);
		}

		void diagonalElementTpem(const TpemSparseType& matrix,
		                         std::vector<RealType> &moment,
		                         size_t ket,
		                         TpemSubspaceType& info,
		                         const RealType& eps) const
		{
			std::vector<RealOrComplexType> tmp(matrix.rank(),0.0);
			std::vector<RealOrComplexType> jm0(matrix.rank(),0.0);
			std::vector<RealOrComplexType> jm1(matrix.rank(),0.0);

			jm0[ket] = 1.0;/* set |j,0> */

			info.clear();
			info.push(ket);

			/* calculate |j,1> = X|j,0> */
			info.sparseProduct(matrix, jm1, jm0,eps);

			RealOrComplexType sum1 = std::conj(jm0[ket]) *  jm1[ket];
			moment[1] += std::real(sum1);

			RealOrComplexType sum2 = 0.0;
			for (size_t p = 0; p < info.top(); p++)
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
			size_t n=moment.size();
			for (size_t m = 2; m < n / 2; m++) {
				/* calculate |tmp> = X|jm1> */
				info.sparseProduct(matrix,tmp, jm1,eps);
				RealOrComplexType sum1 = 0.0;
				RealOrComplexType sum2 = 0.0;
				for (size_t p = 0; p < info.top(); p++) {
					size_t i = info(p);
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

		void  diagonalElementPem(const TpemSparseType& matrix,
		                         std::vector<RealType> &moment,
		                         size_t ket) const
		{
			std::vector<RealOrComplexType> tmp(matrix.rank(),0.0);
			std::vector<RealOrComplexType> jm0(matrix.rank(),0.0);
			std::vector<RealOrComplexType> jm1(matrix.rank(),0.0);

			jm0[ket] = 1.0;	/* set |j,0> */

			/* calculate |j,1> = X|j,0> */
			TpemSubspaceType::sparseProductPem(matrix, jm1, jm0);

			RealOrComplexType sum1 = std::conj(jm0[ket]) * jm1[ket];
			moment[1] += std::real (sum1);

			RealOrComplexType sum2 =  0.0;
			for (size_t j=0;j<matrix.rank();j++)
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
			size_t n=moment.size();
			for (size_t m = 2; m < n / 2; m++) {
				/* calculate |tmp> = X|jm1> */
				TpemSubspaceType::sparseProductPem(matrix, tmp, jm1);
				RealOrComplexType sum1 = 0.0;
				RealOrComplexType sum2 = 0.0;
				for (size_t i=0;i<matrix.rank();i++) {
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

		void subspaceForTrace(TpemSubspaceType& info,
		                      const TpemSparseType& matrix0,
		                      const TpemSparseType& matrix1,
		                      std::vector<RealType>& moment0,
		                      std::vector<RealType>& moment1) const
		{
			const std::vector<size_t>& support = tpemParameters_.support;
			info.clear();

			TpemSubspaceType work(matrix0.rank());
			for (size_t i = 0; i < support.size(); i++) {
				size_t j = support[i];
				diagonalElement(matrix0, moment0, j,work);
				for (size_t p = 0; p < work.top(); p++) info.push(work(p));
				diagonalElement (matrix1, moment1, j,work);
				for (size_t p = 0; p < work.top(); p++) info.push(work(p));
			}
		}

		const TpemParametersType& tpemParameters_; // not the owner, just a ref
		ConcurrencyType& concurrency_;
		CommType comm_;
		GslWrapperType gslWrapper_;
	}; // class Tpem

	template<typename TpemParametersType,typename RealOrComplexType,typename ConcurrencyType>
	const typename ConcurrencyType::CommType
	Tpem<TpemParametersType,RealOrComplexType,ConcurrencyType>::COMM_WORLD = 
	ConcurrencyType::COMM_WORLD; 

} // namespace Tpem

/*@}*/
#endif // TPEM_H
