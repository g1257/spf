
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
#include <string>
#include <vector>
#include "Matrix.h" // in PsimagLite
#include <iostream>

namespace Spf {
	template<typename FieldType_>
	class GeometrySquare45Degrees {
	public:
		//typedef FieldType_ FieldType;
		enum {DIRX=0,DIRY=1,DIRXPY=2,DIRXMY=3};
		
		typedef std::pair<size_t,size_t> PairType;
		
		GeometrySquare45Degrees(size_t l) : l_(l),volume_(2*l*l)
		{
			buildNeighbors();
			std::cout<<(*this);
// 			assert(false);
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

// 		PairType getNeighbour(size_t i,size_t dir) const
// 		{
// 			size_t j = dir*2;
// 			size_t distance = 1;
// 			if (j>=4) {
// 				distance++;
// 				j -= 4;
// 			}
// 			return neighbors_[distance-1](i,j);
// 		}	
// 		
		size_t volume() const { return volume_; }
		
		size_t add(size_t ind,size_t ind2) const
		{
			unimplemented("add");
			return 0;
		}
		
		size_t dim() const { return 2; }
		
		size_t length() const 
		{
			unimplemented("length");
			return 0;
		}

// 		std::string name() const { return "square45degrees"; }

// 		void indexToCoor(std::vector<int> &v,size_t i) const
// 		{
// 			size_t lx = l_;
// 			v[0] = i%lx;
// 			v[1] = size_t(i/lx);
// 		}
		
		size_t coorToIndex(size_t x,size_t y) const
		{
			unimplemented("coorToIndex");
			return 0;
		}

		// direction connecting i and j
		// assume that i and j are n-neighbors or next n-neighbors
		int getDirection(size_t ind,size_t jnd) const
		{
			unimplemented("getDirection");
			return  -1;
		}
// 		
		template<typename T>
		friend std::ostream& operator<<(std::ostream& os,
		                                const GeometrySquare45Degrees<T>& g);

	private:
		
		void unimplemented(const std::string& s) const
		{
			std::string ss = "GeometrySquare45Degrees::" + s;
			ss += " is unimplemented yet (sorry)\n";
			
			throw std::runtime_error(ss.c_str());
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
			size_t n = volume();
			size_t twolx = 2*lx;
			PsimagLite::Matrix<PairType> matrix(n,4);
			for (size_t p=0;p<n;p+=twolx) {
				for (size_t i=p+lx;i<p+twolx;i++) {
					size_t counter = 0;
					int j = i-lx;
					if (j<0) j += n;
					if (i>0 && i%twolx==0) j = i - 1;
					if (i==0) j = n-1;
					matrix(i,counter++) = PairType(j,DIRXMY);
					
					j = i-lx+1;
					if (i+1>=size_t(lx)) {
						size_t k = (i+1);
						if (k%twolx==0) j -= lx;
					}
					if (j<0) j+=n;
					matrix(i,counter++) = PairType(j,DIRXPY);
					
					j = i+lx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRXPY);

					j = i+lx+1;
					if ((i+1)%twolx==0) j -= lx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRXMY);
					
				}
				for (size_t i=p;i<p+lx;i++) {
					size_t counter = 0;
					int j = i-lx-1;
					if (j<0) j += n;
					if (i>0 && i%twolx==0) j = i - 1;
					if (i==0) j = n-1;
					matrix(i,counter++) = PairType(j,DIRXMY);

					j = i-lx;
					if (j<0) j+=n;
					matrix(i,counter++) = PairType(j,DIRXPY);

					j = i+lx-1;
					if (i%twolx==0) j+= lx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRXPY);
					
					j = i+lx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRXMY);
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
			size_t n = volume();
			PsimagLite::Matrix<PairType> matrix(n,4);
			size_t twolx = 2*lx;
			for (size_t p=0;p<n;p+=twolx) {
				for (size_t i=p;i<p+size_t(lx);i++) {
					size_t counter = 0;
					
					int j = i-1;
					if (i%twolx==0) j += lx;

					matrix(i,counter++) = PairType(j,DIRX);

					j = i+1;
					if ((i+1)%lx==0) j -= lx;
					matrix(i,counter++) = PairType(j,DIRX);
					
					j = i - twolx; 
					if (j<0) j += n;
					matrix(i,counter++) = PairType(j,DIRY);
					
					j = i + twolx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRY);
				}
				for (size_t i=p+lx;i<p+twolx;i++) {
					size_t counter = 0;
					int j = i-1;
					size_t k = i-lx;
					if (k%twolx==0) j += lx;
					matrix(i,counter++) = PairType(j,DIRX);

					j = i+1;
					if ((i+1)%lx==0) j -= lx;
					matrix(i,counter++) = PairType(j,DIRX);

					j = i - twolx; 
					if (j<0) j += n;
					matrix(i,counter++) = PairType(j,DIRY);

					j = i + twolx;
					if (size_t(j)>=n) j-=n;
					matrix(i,counter++) = PairType(j,DIRY);
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
		
		size_t g_index(std::vector<int>& x) const
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
		std::vector<PsimagLite::Matrix<PairType> > neighbors_;
	}; //class Square45Degrees
	
	std::ostream& operator<<(std::ostream& os,const std::pair<size_t,size_t>& p)
	{
		os<<p.first<<" "<<p.second;
		return os;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os,const GeometrySquare45Degrees<T>& g)
	{
		for (size_t i=0;i<g.neighbors_.size();i++) {
			os<<"#distanceIndex="<<i<<"\n";
			for (size_t k=0;k<g.neighbors_[i].n_row();k++) {
				os<<"Neighbors of "<<k<<" are ";
				for (size_t l=0;l<g.neighbors_[i].n_col();l++) {
					os<<g.neighbors_[i](k,l).first<<"\t";
				}
				os<<"\n";
			}
		}
		return os;
	}
} // namespace Spf


/*@}*/
#endif // GEOM_SQUARE45DEGREES_H
