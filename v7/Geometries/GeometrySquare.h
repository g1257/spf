
/** \ingroup SPF */
/*@{*/

/*! \file GeometrySquare.h
 *
 * A square lattice
 *
 */
#ifndef GEOM_SQUARE_H
#define GEOM_SQUARE_H
#include "Utils.h"

namespace Spf {
	template<typename FieldType_>
	class GeometrySquare {
		public:
		//typedef FieldType_ FieldType;
		enum {DIRX=0,DIRY=1,DIRXPY=2,DIRXMY=3};
		
		typedef std::pair<size_t,size_t> PairType;
		
		GeometrySquare(size_t l) : l_(l),volume_(l*l)
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
			std::vector<size_t> x(2),y(2);
			index2Coor(x,ind);
			index2Coor(y,ind2);
			for (size_t i=0;i<x.size();i++) {
				x[i] += y[i];
				g_pbc(x[i],l_);
			}
			return g_index(x);
		}
		
		size_t dim() const { return 2; }
		
		private:
		
		void buildNeighbors()
		{
			neighborsAt1();
			neighborsAt2();
		}
		
		void neighborsAt1()
		{
			size_t lx = l_;
			size_t ly = l_;
			size_t zz = 0;
			//PairType zeroVal(0,0);
			psimag::Matrix<PairType> matrix(lx*ly,4);
			for (size_t y=0;y<ly;y++) {
				for (size_t x=0;x<lx;x++) {
					size_t i = x + y*lx;
					size_t xx=x+1; size_t yy=y;
					size_t counter = 0;
					/*if (g_pbc(xx,lvector[0])) {
						border.push_back(0);
					} else {
						border.push_back(-1);
					}*/
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRX);
					xx=x-1;
					/*if (g_pbc(xx,lvector[0])) {
						border.push_back(0);
					} else {
						border.push_back(-1);
					}*/
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRX);
					xx=x; yy=y+1;
					/*if (g_pbc(yy,lvector[1])) {
						border.push_back(1);
					} else {
						border.push_back(-1);
					}*/
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRY);
					yy=y-1;
					/*if (g_pbc(yy,lvector[1])) {
						border.push_back(1);
					} else {
						border.push_back(-1);
					}*/
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRY);
				}
			}
			neighbors_.push_back(matrix);
		}
		
		void neighborsAt2()
		{
			size_t lx = l_;
			size_t ly = l_;
			size_t zz = 0;
			//PairType zeroVal(0,0);
			psimag::Matrix<PairType> matrix(lx*ly,4);
			for (size_t y=0;y<ly;y++) {
				for (size_t x=0;x<lx;x++) {
					size_t i = x + y*lx;
					size_t xx=x+1; size_t yy=y+1;
					size_t counter=0;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXPY);
					xx=x-1; yy=yy-1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXPY);
					xx=x-1; yy=y+1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXMY);
					xx=x+1; yy=y-1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXMY);
				}
			}
			neighbors_.push_back(matrix);
		}
		
		bool g_pbc(size_t& xx, size_t l) const
		{
			int x = xx;
			int L = l;
			bool r=false;
			if (x<0) r=true; 
			if (x>=L) r=true; 
			while(x<0) x+=L;
			while(x>=L) x-=L;
			xx = x;
			return r;
		}
		
		size_t g_index(std::vector<size_t>& x) const
		{
			size_t zz=0;
			return g_index(x[0],x[1],zz);
		}
		
		size_t g_index(size_t& x,size_t& y,size_t& z) const
		{
			size_t lx = l_;
			size_t ly = l_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			//g_pbc(z,lz);
			return x+y*lx; //+z*L*L;
		}

		void index2Coor(std::vector<size_t> &v,size_t i) const
		{
			size_t lx = l_;
			v[0] = i%lx;
			v[1] = size_t(i/lx);
		}
		size_t l_;
		size_t volume_;
		std::vector<psimag::Matrix<PairType> > neighbors_;
	};
	
} // namespace Spf


/*@}*/
#endif
