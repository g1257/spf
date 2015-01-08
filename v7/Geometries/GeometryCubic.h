
/** \ingroup SPF */
/*@{*/

/*! \file GeometryCubic.h
 *
 * A simple cubic lattice
 *
 */
#ifndef GEOM_CUBIC_H
#define GEOM_CUBIC_H

#include "ProgramGlobals.h"
#include "String.h"
#include "Vector.h"

namespace Spf {
	template<typename FieldType_>
	class GeometryCubic {
	public:
		//typedef FieldType_ FieldType;

		static int const DIMENSION = 3;

 		typedef std::pair<SizeType,SizeType> PairType;

 		GeometryCubic(SizeType l) : l_(l),volume_(l*l*l)
		{
			buildNeighbors();
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

		SizeType volume() const { return volume_; }

		SizeType add(SizeType ind,SizeType ind2) const
		{
			PsimagLite::Vector<int>::Type x(2),y(2);
			indexToCoor(x,ind);
			indexToCoor(y,ind2);
			for (SizeType i=0;i<x.size();i++) {
				x[i] += y[i];
				g_pbc(x[i],l_);
			}
			return g_index(x);
		}

		SizeType dim() const { return DIMENSION; }

		SizeType length() const { return l_; }

		template<typename SomeVectorType>
		typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,void>::Type
		indexToCoor(SomeVectorType& v,SizeType i) const
		{
			SizeType lx = l_, ly = l_;
			v[2] = i/(lx*ly);
			SizeType tmp = i - v[2]*lx*ly;
			v[1] = tmp/lx;
			v[0] = tmp % lx;
		}

		template<typename SomeVectorType>
		typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,SizeType>::Type
		coor2Index(const SomeVectorType& v) const
		{
			SomeVectorType length(3,l_);
			int v0 = v[0];
			g_pbc(v0,length[0]);
			SizeType pos=v0;
			SizeType temp=1;

			for (SizeType i=1;i<v.size();i++) {
				int vi = v[i];
				g_pbc(vi,length[i]);
				temp *= length[i-1];
				pos += vi*temp;
			}
			return pos;
		}

		SizeType coorToIndex(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("coorToIndex unimplemented\n");
		}

		bool isBoundary(SizeType i, SizeType j, SizeType d) const
		{
			PsimagLite::Vector<int>::Type vi(3),vj(3);
			SizeType lOver2 = static_cast<SizeType>(l_*0.5);
			indexToCoor(vi,i);
			indexToCoor(vj,j);

			int dist = vi[d] - vj[d];
			SizeType udist = (dist > 0) ? dist : -dist;
			return (udist >= lOver2);
		}

		PsimagLite::String name() const { return "cubic"; }

		SizeType scalarDirection(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("scalarDirection unimplemented\n");
		}

		int getDirection(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("getDirection unimplemented\n");
			return  -1;
		}

		PairType getNeighbour(SizeType,SizeType) const
		{
			throw PsimagLite::RuntimeError("getNeighbour unimplemented\n");
		}

	private:

		void buildNeighbors()
		{
			neighborsAt1();
			//neighborsAt2();
			std::cerr<<
				"WARNING: GeometryCubic: next-nearest neighbors unimplemented\n";
		}

		void neighborsAt1()
		{
			PsimagLite::Matrix<PairType> matrix(volume_,2*DIMENSION);
			for (SizeType i=0;i<volume_;i++) {
				PsimagLite::Vector<SizeType>::Type v(DIMENSION);
				indexToCoor(v,i);
				SizeType x = v[0], y=v[1], z=v[2];

				int zz = z;
				int yy=y;
				SizeType counter = 0;

				int xx=x+1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRX);

				xx=x-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRX);

				xx=x; yy=y+1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRY);

				yy=y-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRY);

				xx = x; yy = y;
				zz = z + 1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRZ);

				zz = z - 1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRZ);
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
			return g_index(x[0],x[1],x[2]);
		}

		SizeType g_index(int& x,int& y,int& z) const
		{
			SizeType lx = l_;
			SizeType ly = l_;
			SizeType lz = l_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			g_pbc(z,lz);
			return x+y*lx+z*lx*ly;
		}

		SizeType l_;
		SizeType volume_;
		PsimagLite::Vector<PsimagLite::Matrix<PairType> >::Type neighbors_;
	}; // class GeometryCubic

} // namespace Spf


/*@}*/
#endif // GEOM_CUBIC_H

