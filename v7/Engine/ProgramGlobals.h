
/** \ingroup SPF */
/*@{*/

/*! \file ProgramGlobals.h

 *
 */
#ifndef SPF_ProgramGlobals_H
#define SPF_ProgramGlobals_H
#include "Map.h"
#include "Vector.h"
#include "PsimagLite.h"

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

	template<typename VectorType, typename VectorType2>
	static
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorType>::True &&
	PsimagLite::IsVectorLike<VectorType2>::True, void>::Type
	vectorPlus(VectorType& a, const VectorType& b, const VectorType2& c)
	{
		const SizeType n = b.size();
		a.resize(n);
		if (n != c.size())
			err("Vector of different sizes\n");
		for (SizeType i = 0; i < n; ++i)
			a[i] = b[i] + c[i];
	}

	template<typename VectorType, typename VectorType2>
	static
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorType>::True &&
	PsimagLite::IsVectorLike<VectorType2>::True, void>::Type
	vectorPlus(VectorType& a, const VectorType& b, const VectorType2& c, const VectorType2& d)
	{
		const SizeType n = b.size();
		if (n != c.size() || n != d.size())
			err("Vector of different sizes\n");
		a.resize(n);
		for (SizeType i = 0; i < n; ++i)
			a[i] = b[i] + c[i] + d[i];
	}

	template<typename VectorType, typename VectorType2>
	static
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorType>::True &&
	PsimagLite::IsVectorLike<VectorType2>::True, void>::Type
	vectorFour(VectorType& a,
	           typename VectorType::value_type bs,
	           const VectorType& b,
	           typename VectorType::value_type cs,
	           const VectorType& c,
	           typename VectorType::value_type ds,
	           const VectorType2& d,
	           typename VectorType::value_type es,
	           const VectorType2& e)
	{
		const SizeType n = b.size();
		if (n != c.size() || n != d.size() || n != e.size())
			err("Vector of different sizes\n");
		a.resize(n);
		for (SizeType i = 0; i < n; ++i)
			a[i] = bs*b[i] + cs*c[i] + ds*d[i] + es*e[i];
	}

	template<typename VectorType, typename VectorType2>
	static
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorType>::True &&
	PsimagLite::IsVectorLike<VectorType2>::True, void>::Type
	vectorMinus(VectorType& a, const VectorType& b, const VectorType2& c)
	{
		const SizeType n = b.size();
		a.resize(n);
		if (n != c.size())
			err("Vector of different sizes\n");
		for (SizeType i = 0; i < n; ++i)
			a[i] = b[i] - c[i];
	}
}; // ProgramGlobals
} // namespace Spf

/*@}*/
#endif // SPF_ProgramGlobals_H
