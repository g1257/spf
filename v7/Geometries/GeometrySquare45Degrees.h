
/** \ingroup SPF */
/*@{*/

/*! \file GeometrySquare45Degrees.h
 *
 * A square lattice
 *
 */
#ifndef GEOM_SQUARE45DEGREES_H
#define GEOM_SQUARE45DEGREES_H
#include <utility>
#include "String.h"
#include <vector>
#include "Matrix.h" // in PsimagLite
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include "ProgramGlobals.h"

namespace Spf {
	template<typename FieldType_>
	class GeometrySquare45Degrees {


		struct SiteInfo {

			enum {TYPE_0,TYPE_1};

			SiteInfo(SizeType indd,SizeType l) : ind(indd)
			{
				init(l);
			}

			SizeType ind;
			SizeType x;
			SizeType y;
			SizeType type;

		private:

			void init(SizeType l)
			{
				SizeType lx = l;
				SizeType twolx = 2*lx;
				div_t divresult = div(ind,twolx);

				if (SizeType(divresult.rem)<lx) {
					type = TYPE_0;
					x = divresult.rem;
					y = divresult.quot;
					return;
				}
				type = TYPE_1;
				x = divresult.rem-lx;
				y = divresult.quot;
			}
		};

		typedef SiteInfo SiteInfoType;

	public:

		typedef std::pair<SizeType,SizeType> PairType;

		GeometrySquare45Degrees(SizeType l) : l_(l),volume_(2*l*l)
		{
			buildNeighbors();
			std::cout<<(*this);
// 			assert(false);
		}

		SizeType distances() const { return neighbors_.size(); }

		SizeType z(SizeType distance=1) const
		{
			return neighbors_[distance-1].n_col();
		}

		// j-th neighbor of i at distance (starts from 1 for compatibility)
		PairType neighbor(SizeType i,SizeType j,SizeType distance=1) const
		{
			return neighbors_[distance-1](i,j);
		}

// 		PairType getNeighbour(SizeType i,SizeType dir) const
// 		{
// 			SizeType j = dir*2;
// 			SizeType distance = 1;
// 			if (j>=4) {
// 				distance++;
// 				j -= 4;
// 			}
// 			return neighbors_[distance-1](i,j);
// 		}
//
		SizeType volume() const { return volume_; }

		SizeType add(SizeType ind,SizeType ind2) const
		{
			SiteInfoType s1(ind,l_);
			SiteInfoType s2(ind2,l_);

			SizeType lx = l_;
			SizeType ly = l_;
			SizeType xsum = s1.x + s2.x;
			SizeType ysum = s1.y + s2.y;
			if (xsum >= lx) xsum -= lx;
			if (ysum >= ly) ysum -= ly;

			if (s1.type==s2.type) {
				SizeType shift = (s1.type==SiteInfoType::TYPE_0) ? 0 : 1;
				return siteAt(xsum+shift,ysum+shift);
			}
			SizeType j = siteAt(xsum,ysum) + lx;
			SizeType n = volume();
			if (SizeType(j)>=n) j-=n;
			return j;
		}

		SizeType dim() const { return 2; }

		SizeType length() const
		{
			unimplemented("length");
			return 0;
		}

 		PsimagLite::String name() const { return "square45degrees"; }

// 		void indexToCoor(PsimagLite::Vector<int>::Type& v,SizeType i) const
// 		{
// 			SizeType lx = l_;
// 			v[0] = i%lx;
// 			v[1] = SizeType(i/lx);
// 		}

		SizeType coorToIndex(SizeType,SizeType) const
		{
			unimplemented("coorToIndex");
			return 0;
		}

		// direction connecting i and j
		// assume that i and j are n-neighbors or next n-neighbors
		int getDirection(SizeType,SizeType) const
		{
			unimplemented("getDirection");
			return  -1;
		}

		SizeType scalarDirection(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("scalarDirection unimplemented\n");
		}

		template<typename SomeVectorType>
		typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,SizeType>::Type
		coor2Index(const SomeVectorType& v) const
		{
			throw PsimagLite::RuntimeError("coor2Index unimplemented\n");
		}

		PairType getNeighbour(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("getNeighbour unimplemented\n");
		}

		template<typename SomeVectorType>
		typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,void>::Type
		indexToCoor(SomeVectorType&,SizeType) const
		{
			throw PsimagLite::RuntimeError("indexToCoor unimplemented\n");
		}

		bool isBoundary(SizeType, SizeType, SizeType) const
		{
			throw PsimagLite::RuntimeError("isBoundary unimplemented\n");
		}

		template<typename T>
		friend std::ostream& operator<<(std::ostream& os,
		                                const GeometrySquare45Degrees<T>& g);

	private:



		SizeType siteAt(SizeType x,SizeType y) const
		{
			SizeType lx = l_;
			SizeType ly = l_;
			SizeType twolx = 2*lx;
			if (x>=lx) x-=lx;
			if (y>=ly) y-=ly;
			return y*twolx + x;
		}

		void unimplemented(const PsimagLite::String& s) const
		{
			PsimagLite::String ss = "GeometrySquare45Degrees::" + s;
			ss += " is unimplemented yet (sorry)\n";

			throw PsimagLite::RuntimeError(ss.c_str());
		}

		bool periodicEqualTo(int a,int b) const
		{
			if (a==b) return true;
			int c = -b*(l_-1);
			if (a==c) return true;
			return false;
		}

		void buildNeighbors()
		{
			neighborsAt1();
			neighborsAt2();
		}

		void neighborsAt1()
		{
			int lx = l_;
// 			int ly = l_;
// 			int zz = 0;
			//PairType zeroVal(0,0);
			SizeType n = volume();
			SizeType twolx = 2*lx;
			PsimagLite::Matrix<PairType> matrix(n,4);
			for (SizeType p=0;p<n;p+=twolx) {
				for (SizeType i=p+lx;i<p+twolx;i++) {
					SizeType counter = 0;
					int j = i-lx;
					if (j<0) j += n;
					if (i>0 && i%twolx==0) j = i - 1;
					if (i==0) j = n-1;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXMY);

					j = i-lx+1;
					if (i+1>=SizeType(lx)) {
						SizeType k = (i+1);
						if (k%twolx==0) j -= lx;
					}
					if (j<0) j+=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXPY);

					j = i+lx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXPY);

					j = i+lx+1;
					if ((i+1)%twolx==0) j -= lx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXMY);

				}
				for (SizeType i=p;i<p+lx;i++) {
					SizeType counter = 0;
					int j = i-lx-1;
					if (j<0) j += n;
					if (i>0 && i%twolx==0) j = i - 1;
					if (i==0) j = n-1;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXMY);

					j = i-lx;
					if (j<0) j+=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXPY);

					j = i+lx-1;
					if (i%twolx==0) j+= lx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXPY);

					j = i+lx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRXMY);
				}
			}
			neighbors_.push_back(matrix);
		}

		void neighborsAt2()
		{
			int lx = l_;
// 			int ly = l_;
// 			int zz = 0;
			//PairType zeroVal(0,0);
			SizeType n = volume();
			PsimagLite::Matrix<PairType> matrix(n,4);
			SizeType twolx = 2*lx;
			for (SizeType p=0;p<n;p+=twolx) {
				for (SizeType i=p;i<p+SizeType(lx);i++) {
					SizeType counter = 0;

					int j = i-1;
					if (i%twolx==0) j += lx;

					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRX);

					j = i+1;
					if ((i+1)%lx==0) j -= lx;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRX);

					j = i - twolx;
					if (j<0) j += n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRY);

					j = i + twolx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRY);
				}
				for (SizeType i=p+lx;i<p+twolx;i++) {
					SizeType counter = 0;
					int j = i-1;
					SizeType k = i-lx;
					if (k%twolx==0) j += lx;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRX);

					j = i+1;
					if ((i+1)%lx==0) j -= lx;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRX);

					j = i - twolx;
					if (j<0) j += n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRY);

					j = i + twolx;
					if (SizeType(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,ProgramGlobals::DIRY);
				}
			}
			neighbors_.push_back(matrix);
		}

		bool g_pbc(int& x, SizeType l) const
		{
			int L = l;
			bool r=false;
			if (x<0) r=true;
			if (x>=L) r=true;
			while(x<0) x+=L;
			while(x>=L) x-=L;
			return r;
		}

		SizeType g_index(PsimagLite::Vector<int>::Type& x) const
		{
			int zz=0;
			return g_index(x[0],x[1],zz);
		}

		SizeType g_index(int& x,int& y,int& z) const
		{
			SizeType lx = l_;
			SizeType ly = l_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			//g_pbc(z,lz);
			return x+y*lx; //+z*L*L;
		}

		SizeType l_;
		SizeType volume_;
		PsimagLite::Vector<PsimagLite::Matrix<PairType> >::Type neighbors_;
	}; //class Square45Degrees

	template<typename T>
	std::ostream& operator<<(std::ostream& os,const GeometrySquare45Degrees<T>& g)
	{
		for (SizeType i=0;i<g.neighbors_.size();i++) {
			os<<"#distanceIndex="<<i<<"\n";
			for (SizeType k=0;k<g.neighbors_[i].n_row();k++) {
				os<<"Neighbors of "<<k<<" are ";
				for (SizeType l=0;l<g.neighbors_[i].n_col();l++) {
					os<<g.neighbors_[i](k,l).first<<"\t";
				}
				os<<"\n";
			}
		}
		os<<"-------------------------\n";
		for (SizeType i=0;i<g.volume();i++) {
			for (SizeType j=0;j<g.volume();j++) {
				SizeType k = g.add(i,j);
				os<<i<<"\t+\t"<<j<<"\t=\t"<<k<<"\n";
			}
		}
		return os;
	}
} // namespace Spf


/*@}*/
#endif // GEOM_SQUARE45DEGREES_H
