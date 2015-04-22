
/** \ingroup SPF */
/*@{*/

/*! \file FakeParams.h
 *
 *
 *
 */
#ifndef FAKE_PARAMS_H
#define FAKE_PARAMS_H
#include "AllocatorCpu.h"

namespace Spf {

struct FakeParams {

	FakeParams(PsimagLite::String dynvarsfile1,
	           int long long randomSeed1,
	           PsimagLite::String options1)
	    : dynvarsfile(dynvarsfile1),randomSeed(randomSeed1),options(options1)
	{}

	PsimagLite::String dynvarsfile;
	int long long randomSeed;
	PsimagLite::String options;
}; // FakeParams
} // namespace Spf

/*@}*/
#endif // FAKE_PARAMS_H
