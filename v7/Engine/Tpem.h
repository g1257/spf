
/** \ingroup SPF */
/*@{*/

/*! \file Tpem.h
 *
 *  
 *
 */
#ifndef TPEM_H
#define TPEM_H
#include <gsl/gsl_integration.h>
#include "CrsMatrix.h"

namespace Spf {

	template<typename BaseFunctorType>
	class Tpem {
		typedef typename BaseFunctorType::RealType RealType;
	public:
		
		typedef PsimagLite::CrsMatrix<RealType> TpemSparseType;

		struct MyFunctionParams {
			MyFunctionParams(const BaseFunctorType& functor1)
			: functor(functor1) { }

			const BaseFunctorType& functor;
			size_t m;
		};

		typedef MyFunctionParams MyFunctionParamsType;

		Tpem()
		{}

		void calcCoeffs(std::vector<RealType> &vobs,
		                const BaseFunctorType& obsFunc)
		{
			std::vector<RealType> pts(2);
			pts[0]= -1.0;
			pts[1] = 1.0;
			int npts = pts.size();
			RealType epsabs=1e-9;
			RealType epsrel=1e-9;
			
			int limit = 1e6;
			gsl_integration_workspace *workspace= gsl_integration_workspace_alloc(limit+2);
			
			RealType result = 0,abserr = 0;
			
			gsl_function f;
			f.function= &Tpem<BaseFunctorType>::myFunction;
			MyFunctionParamsType params(obsFunc);
			f.params = &params;
			
			for (size_t m=0;m<cutoff_;m++) {
				params.m = m;
				gsl_integration_qagp(&f,&(pts[0]),npts,epsabs,epsrel,limit,workspace,&result,&abserr);
				vobs[m] = result;
			}
			gsl_integration_workspace_free (workspace);
		}
		
		double tpem_expansion (vector<double> const & moments, vector<double> const &coeffs) 
		{
			unsigned int	i,n=moments.size();
			double	ret = 0.0;
			
			if (n!=coeffs.size()) {
				cerr<<"tpem_expansion: coeffs and moments of different size\n";
			}
			for (i = 0; i < n; i++)
				ret += moments[i] * coeffs[i];
			return ret;
		}
	private:
		
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

		static RealType chebyshev(int m,RealType x) {
			RealType tmp;
			int p;
			if (m==0) return 1;
			if (m==1) return x;
			
			if ((m%2)==0) {
				p=m/2;
				tmp=chebyshev(p,x);
				return (2*tmp*tmp-1);
			}
			else {
				p=(m-1)/2;
				return (2*tpem_chebyshev(p,x)*tpem_chebyshev(p+1,x)-x);
			}
		}
	}; // class Tpem
} // namespace Spf

/*@}*/
#endif // TPEM_H
