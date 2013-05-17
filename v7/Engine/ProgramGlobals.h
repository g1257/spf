
/** \ingroup SPF */
/*@{*/

/*! \file ProgramGlobals.h

 *
 */
#ifndef SPF_ProgramGlobals_H
#define SPF_ProgramGlobals_H
#include "String.h"
#include "Map.h"

namespace Spf {

class ProgramGlobals {

public:

	enum {DIRX=0,DIRY=1,DIRXPY=2,DIRXMY=3,DIRZ=0};

	template<typename MapType>
	static typename PsimagLite::HasType<PsimagLite::IsMapLike<MapType>::True,void>::Type
	checkMcWindow(const MapType& mcWindow,const PsimagLite::String& arg)
	{
		if (mcWindow.find(arg.c_str()) != mcWindow.end()) return;

			PsimagLite::String str("checkMcWindow(): MonteCarloWindow[");
			str += arg + "] needed in input file\n";

			throw PsimagLite::RuntimeError(str);
	}

	template<typename T>
	static T square(const T& t)
	{
		return t*t;
	}
}; // ProgramGlobals
} // namespace Spf

/*@}*/
#endif // SPF_ProgramGlobals_H
