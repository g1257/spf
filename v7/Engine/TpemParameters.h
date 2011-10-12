
/** \ingroup SPF */
/*@{*/

/*! \file TpemParameters.h
 *
 *  
 *
 */
#ifndef TPEM_PARAMETERS_H
#define TPEM_PARAMETERS_H
#include <gsl/gsl_integration.h>
#include "CrsMatrix.h"
#include <cassert>

namespace Spf {
	struct TpemParameters {
		enum {TPEM,PEM};
		size_t cutoff;
		size_t algorithm;
	}; // struct TpemParameters
} // namespace Spf

/*@}*/
#endif // TPEM_PARAMETERS_H
