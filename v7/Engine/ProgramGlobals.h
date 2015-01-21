
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

	template<typename SomeVectorType>
	static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	void>::Type checkMcWindow(const SomeVectorType& mcWindow, SizeType i)
	{
		if (i < mcWindow.size()) return;

			PsimagLite::String str("checkMcWindow(): MonteCarloWindow must have");
			str += "at least " + ttos(i+1) + " entries.\n";

			throw PsimagLite::RuntimeError(str);
	}

	template<typename T>
	static void printGeometry(std::ostream& os,const T& g)
	{
		SizeType n = g.volume();

		for (SizeType d=0;d<g.distances();++d) {
			os<<"#Distance Number="<<d<<"\n";
			for (SizeType i=0;i<n;++i) {
				std::cout<<"Neighbors of "<<i<<" are ";
				for (SizeType k=0;k<g.z(d+1);k++) {
					SizeType j = g.neighbor(i,k,d+1).first;
					os<<j<<" ";
				}

				os<<"\n";
			}
		}
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
