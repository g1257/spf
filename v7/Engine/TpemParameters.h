
/** \ingroup TPEM */
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

namespace Tpem {
	template<typename IoInType,typename MuBetaStructType>
	struct TpemParameters {
		typedef typename MuBetaStructType::RealType RealType;
		enum {TPEM,PEM};

		TpemParameters(IoInType& io,const MuBetaStructType& engineParams)
		: mu(engineParams.mu),beta(engineParams.beta)
		{
			io.readline(algorithm,"TpemAlgorithm=");
			io.readline(cutoff,"TpemCutoff=");
			io.readline(eps,"TpemEps=");
			io.readline(a,"TpemA=");
			io.readline(b,"TpemB=");
			io.read(support,"TpemSupport");
		}

		const RealType& mu;
		const RealType& beta;
		size_t cutoff;
		size_t algorithm;
		RealType eps;
		RealType a;
		RealType b;
		std::vector<size_t> support;
	}; // struct TpemParameters
} // namespace Tpem

/*@}*/
#endif // TPEM_PARAMETERS_H
