
/** \ingroup SPF */
/*@{*/

/*! \file GeometryCubic.h
 *
 * A simple cubic lattice
 *
 */
#ifndef GEOM_CUBIC_H
#define GEOM_CUBIC_H
#include "Utils.h"

namespace Spf {
	template<typename FieldType_>
	class GeometryCubic {
	public:
		//typedef FieldType_ FieldType;
		enum {DIRX=0,DIRY=1,DIRZ=0};
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
			std::vector<int> x(2),y(2);
			index2Coor(x,ind);
			index2Coor(y,ind2);
			for (size_t i=0;i<x.size();i++) {
				x[i] += y[i];
				g_pbc(x[i],l_);
			}
			return g_index(x);
		}
		
		size_t dim() const { return DIMENSION; }
		
		size_t length() const { return l_; }
		
		void index2Coor(std::vector<size_t> &v,size_t i) const
		{
			size_t lx = l_, ly = l_;
			v[2] = i/(lx*ly);
			size_t tmp = i - v[2]*lx*ly;
			v[1] = tmp/lx;
			v[0] = tmp % lx;
		}

		std::string name() const { return "cubic"; }

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
				std::vector<size_t> v(DIMENSION);
				index2Coor(v,i);
				size_t x = v[0], y=v[1], z=v[2];

				int zz = z;
				int yy=y;
				size_t counter = 0;

				int xx=x+1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRX);

				xx=x-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRX);

				xx=x; yy=y+1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRY);

				yy=y-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRY);

				xx = x; yy = y;
				zz = z + 1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRZ);

				zz = z - 1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRZ);
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
		
		size_t g_index(std::vector<int>& x) const
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
		std::vector<PsimagLite::Matrix<PairType> > neighbors_;
	}; // class GeometryCubic
	
} // namespace Spf


/*@}*/
#endif // GEOM_CUBIC_H

