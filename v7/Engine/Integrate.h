
/** \ingroup SPF */
/*@{*/

/*!
 *
 *
 *
 */
#ifndef INTEGRATE_SPF_H
#define INTEGRATE_SPF_H
#include "GslWrapper.h"
#include "Vector.h"

namespace Spf {

template<typename FunctionParamsType>
double myFunction(double x, void * p)
{
	FunctionParamsType* params = (FunctionParamsType*)p;

	return params->operator()(x);
}

template<typename SomeFunctorType>
typename SomeFunctorType::RealType integrate(SomeFunctorType& functor1,
                                             typename SomeFunctorType::RealType x1,
                                             typename SomeFunctorType::RealType x2)
{
	typedef PsimagLite::GslWrapper GslWrapperType;
	typedef typename SomeFunctorType::RealType RealType;
	typedef typename SomeFunctorType::ParametersType ParametersType;

	typename PsimagLite::Vector<RealType>::Type pts(2);
	pts[0] = x1;
	pts[1] = x2;
	RealType epsabs=1e-9;
	RealType epsrel=1e-9;

	GslWrapperType gslWrapper;

	size_t limit = 1000000;
	GslWrapperType::gsl_integration_workspace *workspace =
	        gslWrapper.gsl_integration_workspace_alloc(limit+2);

	RealType result = 0;
	RealType abserr = 0;

	GslWrapperType::gsl_function f;
	f.params = &functor1;
	f.function= myFunction<SomeFunctorType>;
	gslWrapper.gsl_integration_qagp(&f,&(pts[0]),pts.size(),epsabs,epsrel,limit,
	                                 workspace,&result,&abserr);

	gslWrapper.gsl_integration_workspace_free (workspace);
	return result;
}

} // namespace Spf

/*@}*/
#endif // INTEGRATE_SPF_H
