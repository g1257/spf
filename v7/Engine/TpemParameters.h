
/** \ingroup SPF */
/*@{*/

/*! \file TpemParameters.h
 *
 *  Do not add functions here except a ctor
 *
 */
#ifndef TPEM_PARAMETERS_H
#define TPEM_PARAMETERS_H
#include <gsl/gsl_integration.h>
#include "CrsMatrix.h"
#include <cassert>

namespace Spf {
	template<typename IoInType,typename RealType_>
	struct TpemParameters {
		typedef RealType_ RealType;
		enum {TPEM,PEM};

		TpemParameters(IoInType& io)
		{
			io.readline(algorithm,"TpemAlgorithm=");
			io.readline(cutoff,"TpemCutoff=");
			io.readline(a,"TpemA=");
			io.readline(b,"TpemB=");
		}

		size_t cutoff;
		size_t algorithm;
		RealType a;
		RealType b;
	}; // struct TpemParameters
} // namespace Spf

/*@}*/
#endif // TPEM_PARAMETERS_H
