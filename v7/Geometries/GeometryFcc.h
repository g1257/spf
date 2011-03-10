
/** \ingroup SPF */
/*@{*/

/*! \file GeometryFcc.h
 *
 * A FCC lattice
 *
 */
#ifndef GEOM_FCC_H
#define GEOM_FCC_H
#include "GeometryCubic.h"
#include "Vector.h"
#include "TypeToString.h"

namespace Spf {
	template<typename RealType>
	class GeometryFcc {
	public:
		enum {DIRX=0,DIRY=1,DIRZ=0};
		static size_t const DIMENSION = 3;
		static size_t const COORDINATION = 12;
		static size_t const BASIS_FOR_CUBIC  = 4;
		const RealType MIN_DISTANCE;

 		typedef std::pair<size_t,size_t> PairType;
 		typedef GeometryCubic<RealType> GeometryCubicType;

 		GeometryFcc(size_t l) : MIN_DISTANCE(0.5) , l_(l),
 				cube_(l), volume_(cube_.volume()*BASIS_FOR_CUBIC)
		{
			buildNeighbors();
			throw std::runtime_error("FCC DIRECTIONS OF NEIGHBORS!\n");
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
			throw std::runtime_error("FCC ADD\n");
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
		
	private:
		
		void buildNeighbors()
		{
			initBasis();
			neighborsAt1();
		}
		
		void neighborsAt1()
		{
			PsimagLite::Matrix<PairType> matrix(volume_,COORDINATION);
			for (size_t i=0;i<volume_;i++) {
				//! find sites that are a distance min from i
				//! and store them on tmpVec
				std::vector<size_t> tmpVec;
				findNn(tmpVec,i);
				size_t counter = 0;
				for (size_t j=0;j<tmpVec.size();j++)
					matrix(i,counter++) = PairType(tmpVec[j],DIRX);
			}
			neighbors_.push_back(matrix);
		}
		
		void findNn(std::vector<size_t> &v,size_t site) const
		{
			for (size_t i=0;i<volume_;i++) {
				RealType dd=dist(i,site);
				if (fabs(dd-MIN_DISTANCE)<1e-6) {
					v.push_back(i);
				}
			}
			if (v.size()!=COORDINATION) {
				std::string s ="GeometryFcc: "+PsimagLite::typeToString(site) +
					" has " + PsimagLite::typeToString(v.size()) +
					" neighbours but it should be " +
					PsimagLite::typeToString(COORDINATION) +
						" instead\n";
				throw std::runtime_error(s.c_str());
			}
		}

		RealType dist(size_t site1, size_t site2) const
		{
			std::vector<RealType> r1,r2;

			index2Coor(r1,site1);
			index2Coor(r2,site2);
			return distance(r1,r2);
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

		//! order of the 4 basis vectors is
		//! (0,0,0) (0,a/2,a/2), (a/2,0,a/2), (a/2,a/2,0)
		//! coordinates are (x,y,z,b) where b is the basis
		void index2Coor(std::vector<RealType> &r,size_t site) const
		{
			size_t b = site/cube_.volume();
			size_t i = site % cube_.volume();
			std::vector<size_t> rcube(DIMENSION);
			cube_.index2Coor(rcube,i);
			r = basisVector_[b] + rcube;
		}

		RealType distance(
				const std::vector<RealType>& r1,
				const std::vector<RealType>& r2) const
		{
			std::vector<RealType> r(r1.size());

			for (size_t i=0;i<r1.size();i++) {
				r[i]=r1[i]-r2[i];
				while(r[i]<MIN_DISTANCE) r[i]+=l_;
				while(r[i]>=l_-MIN_DISTANCE) r[i]-=l_;
			}
			return r*r;
		}

		void initBasis()
		{
			std::vector<RealType> tmpVector(DIMENSION,0);

			basisVector_.push_back(tmpVector);

			tmpVector[0]=0; tmpVector[1]=0.5; tmpVector[2]=0.5;
			basisVector_.push_back(tmpVector);

			tmpVector[0]=0.5; tmpVector[1]=0; tmpVector[2]=0.5;
			basisVector_.push_back(tmpVector);

			tmpVector[0]=0.5; tmpVector[1]=0.5; tmpVector[2]=0;
			basisVector_.push_back(tmpVector);
		}

		size_t l_;
		GeometryCubicType cube_;
		size_t volume_;
		std::vector<std::vector<RealType> > basisVector_;
		std::vector<PsimagLite::Matrix<PairType> > neighbors_;
	}; // class GeometryFcc
	
} // namespace Spf


/*@}*/
#endif // GEOM_FCC_H

