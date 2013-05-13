
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

namespace Spf {
	template<typename FieldType_>
	class GeometryCubic {
	public:
		//typedef FieldType_ FieldType;

		static int const DIMENSION = 3;
		
 		typedef std::pair<size_t,size_t> PairType;
		
 		GeometryCubic(size_t l) : l_(l),volume_(l*l*l)
		{
			buildNeighbors();
		}
		
		size_t z(size_t distance=1) const
		{
			return neighbors_[distance-1].n_col();
		}
		
		// j-th neighbor of i at distance (starts from 1 for compatibility)
		PairType neighbor(size_t i,size_t j,size_t distance=1) const
		{
			return neighbors_[distance-1](i,j);
		}
		
		size_t volume() const { return volume_; }
		
		size_t add(size_t ind,size_t ind2) const
		{
			PsimagLite::Vector<int>::Type x(2),y(2);
			indexToCoor(x,ind);
			indexToCoor(y,ind2);
			for (size_t i=0;i<x.size();i++) {
				x[i] += y[i];
				g_pbc(x[i],l_);
			}
			return g_index(x);
		}
		
		size_t dim() const { return DIMENSION; }
		
		size_t length() const { return l_; }
		
		void indexToCoor(PsimagLite::Vector<size_t>::Type& v,size_t i) const
		{
			size_t lx = l_, ly = l_;
			v[2] = i/(lx*ly);
			size_t tmp = i - v[2]*lx*ly;
			v[1] = tmp/lx;
			v[0] = tmp % lx;
		}

		size_t coor2Index(const PsimagLite::Vector<size_t>::Type& v) const
		{
			PsimagLite::Vector<size_t>::Type length(3,l_);
			int v0 = v[0];
			g_pbc(v0,length[0]);
			size_t pos=v0;
			size_t temp=1;

			for (size_t i=1;i<v.size();i++) {
				int vi = v[i];
				g_pbc(vi,length[i]);
				temp *= length[i-1];
				pos += vi*temp;
			}
			return pos;
		}

		PsimagLite::String name() const { return "cubic"; }

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
			for (size_t i=0;i<volume_;i++) {
				PsimagLite::Vector<size_t>::Type v(DIMENSION);
				indexToCoor(v,i);
				size_t x = v[0], y=v[1], z=v[2];

				int zz = z;
				int yy=y;
				size_t counter = 0;

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
		
		bool g_pbc(int& x, size_t l) const
		{
			int L = l;
			bool r=false;
			if (x<0) r=true; 
			if (x>=L) r=true; 
			while(x<0) x+=L;
			while(x>=L) x-=L;
			return r;
		}
		
		size_t g_index(PsimagLite::Vector<int>::Type& x) const
		{
			return g_index(x[0],x[1],x[2]);
		}
		
		size_t g_index(int& x,int& y,int& z) const
		{
			size_t lx = l_;
			size_t ly = l_;
			size_t lz = l_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			g_pbc(z,lz);
			return x+y*lx+z*lx*ly;
		}

		size_t l_;
		size_t volume_;
		PsimagLite::Vector<PsimagLite::Matrix<PairType> >::Type neighbors_;
	}; // class GeometryCubic
	
} // namespace Spf


/*@}*/
#endif // GEOM_CUBIC_H

