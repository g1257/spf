
/** \ingroup SPF */
/*@{*/

/*! \file GeometrySquare.h
 *
 * A square lattice
 *
 */
#ifndef GEOM_SQUARE_H
#define GEOM_SQUARE_H
#include <utility>
#include "String.h"
#include <vector>
#include "Matrix.h" // in PsimagLite

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

		PairType getNeighbour(size_t i,size_t dir) const
		{
			size_t j = dir*2;
			size_t distance = 1;
			if (j>=4) {
				distance++;
				j -= 4;
			}
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
		
		size_t dim() const { return 2; }
		
		size_t length() const { return l_; }
		
		PsimagLite::String name() const { return "square"; }

		void indexToCoor(PsimagLite::Vector<int>::Type& v,size_t i) const
		{
			size_t lx = l_;
			v[0] = i%lx;
			v[1] = size_t(i/lx);
		}
		
		size_t coorToIndex(size_t x,size_t y) const
		{
			size_t lx = l_;
			return x + y*lx;
		}

		// direction connecting i and j
		// assume that i and j are n-neighbors or next n-neighbors
		int getDirection(size_t ind,size_t jnd) const
		{
			PsimagLite::Vector<int>::Type vi(2),vj(2);
			indexToCoor(vi,ind);
			indexToCoor(vj,jnd);
			
			for (size_t i=0;i<vi.size();i++) {
				vi[i] -= vj[i];
			}
			
			if (vi[0] == 0 && periodicEqualTo(vi[1],1)) return DIRY;
			if (vi[0] == 0 && periodicEqualTo(vi[1],-1)) return DIRY;
			
			if (vi[1] == 0 &&  periodicEqualTo(vi[0],1)) return DIRX;
			if (vi[1] == 0 &&  periodicEqualTo(vi[0],-1)) return DIRX;
			
			if (periodicEqualTo(vi[0],1) && periodicEqualTo(vi[1],1)) return DIRXPY;
			if (periodicEqualTo(vi[0],-1) && periodicEqualTo(vi[1],-1)) return DIRXPY;
			
			if (periodicEqualTo(vi[0],1) && periodicEqualTo(vi[1],-1)) return DIRXMY;
			if (periodicEqualTo(vi[0],-1) && periodicEqualTo(vi[1],1)) return DIRXMY;
			
			return -1;
		}
		
		template<typename T>
		friend std::ostream& operator<<(std::ostream& os,
		                                const GeometrySquare<T>& g);

	private:
		
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
			int ly = l_;
			int zz = 0;
			//PairType zeroVal(0,0);
			PsimagLite::Matrix<PairType> matrix(lx*ly,4);
			for (int y=0;y<ly;y++) {
				for (int x=0;x<lx;x++) {
					size_t i = x + y*lx;
					int xx=x+1;
					int yy=y;
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
			int lx = l_;
			int ly = l_;
			int zz = 0;
			//PairType zeroVal(0,0);
			PsimagLite::Matrix<PairType> matrix(lx*ly,4);
			for (int y=0;y<ly;y++) {
				for (int x=0;x<lx;x++) {
					size_t i = x + y*lx;
					int xx=x+1; 
					int yy=y+1;
					size_t counter=0;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXPY);
					xx=x-1; yy=y-1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXPY);
					xx=x-1; yy=y+1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXMY);
					xx=x+1; yy=y-1;
					matrix(i,counter++) = PairType(g_index(xx,yy,zz),DIRXMY);
				}
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
			int zz=0;
			return g_index(x[0],x[1],zz);
		}

		size_t g_index(int& x,int& y,int& z) const
		{
			size_t lx = l_;
			size_t ly = l_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			//g_pbc(z,lz);
			return x+y*lx; //+z*L*L;
		}

		size_t l_;
		size_t volume_;
		PsimagLite::Vector<PsimagLite::Matrix<PairType> >::Type neighbors_;
	};
	
	std::ostream& operator<<(std::ostream& os,const std::pair<size_t,size_t>& p)
	{
		os<<p.first<<" "<<p.second;
		return os;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os,const GeometrySquare<T>& g)
	{
		for (size_t i=0;i<g.neighbors_.size();i++) {
			os<<"#i="<<i<<"\n";
			for (size_t k=0;k<g.neighbors_[i].n_row();k++) {
				for (size_t l=0;l<g.neighbors_[i].n_col();l++) {
					os<<g.neighbors_[i](k,l).first<<" ";
				}
				os<<"\n";
			}
		}
		return os;
	}
} // namespace Spf


/*@}*/
#endif
