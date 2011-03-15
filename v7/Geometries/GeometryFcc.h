
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
		
		//! coordinates are (x,y,z,b) where b is the basis
		size_t add(size_t ind,size_t ind2) const
		{
			// transform ind --> rvector, bvector and ind2 --> r2vector, b2vector
			size_t b1 = ind/volume_;
			size_t indr = ind % volume_;
			std::vector<size_t> rvector(3);
			cube_.index2Coor(rvector,indr);

			size_t b2 = ind2/volume_;
			size_t ind2r = ind2 % volume_;
			std::vector<size_t> r2vector(3);
			cube_.index2Coor(r2vector,ind2r);

			// sum these two and obtain rsumvector bsumvector
			std::vector<RealType> rsumvector(3),bsumvector(3);
			add(rsumvector,bsumvector,rvector,b1,r2vector,b2);

			// transform rsumvector and bsumvector into an index that is returned
			int b = utils::isInVector(basisVector_,bsumvector);
			if (b<0) {
				std::string s = "Fcc::add(...) INTERNAL ERROR.\n";
				s += "Exiting at this point " + std::string(__FILE__) +
				" "+ ttos(__LINE__) + "\n";
				throw std::runtime_error(s.c_str());
			}
			std::vector<size_t> r(3);
			for (size_t i=0;i<rsumvector.size();i++) r[i] = rsumvector[i];
			indr = cube_.coor2Index(r);
			return indr + size_t(b)*cube_.volume();
		}
		
		size_t scalarDirection(size_t site1,size_t site2) const
		{
			std::vector<RealType> r(3);
			direction(r,site1,site2);
			std::vector<RealType> r2(3);
			scDirAux(r2,r);
			return size_t(r2[0]+2*r2[1]+4*r2[2]);
		}

		size_t dim() const { return DIMENSION; }
		
		size_t length() const { return l_; }
		
		std::string name() const { return "fcc"; }

	private:
		
		void buildNeighbors()
		{
			initBasis();
			neighborsAt1();
		}
		
		void neighborsAt1()
		{
			PsimagLite::Matrix<PairType> matrix(volume_,COORDINATION);
			std::vector<RealType> tmpv(DIMENSION);
			std::vector<RealType> rvector(DIMENSION);

			for (size_t i=0;i<volume_;i++) {
				//! find sites that are a distance min from i
				//! and store them on tmpVec
				std::vector<size_t> tmpVec;
				findNn(tmpVec,i);
				size_t counter = 0;
				for (size_t j=0;j<tmpVec.size();j++) {
					direction(rvector,i,tmpVec[j]);
					scDirAux(tmpv,rvector);
					size_t dir = size_t(tmpv[0]+2*tmpv[1]+4*tmpv[2]);
					matrix(i,counter++) = PairType(tmpVec[j],dir);
				}
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
			modifiedDifference(r,r1,r2);
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

		void direction(
				std::vector<RealType> &r,
				size_t site1,
				size_t site2) const
		{
			std::vector<RealType> r1,r2;

			index2Coor(r1,site1);
			index2Coor(r2,site2);

			modifiedDifference(r,r1,r2);
		}

		void modifiedDifference(
				std::vector<RealType>& r,
				const std::vector<RealType>& r1,
				const std::vector<RealType>& r2) const
		{

			for (size_t i=0;i<r1.size();i++) {
				r[i]=r1[i]-r2[i];
				while(r[i]<MIN_DISTANCE) r[i]+=l_;
				while(r[i]>=l_-MIN_DISTANCE) r[i]-=l_;
			}
		}

		void scDirAux(
				std::vector<RealType>& v,
				const std::vector<RealType>& r) const
		{
			if (r[0]==0) {
				v[2]=0;
				v[0]=r[1]+MIN_DISTANCE;
				v[1]=r[2]+MIN_DISTANCE;
			} else if (r[1]==0) {
				v[2]=1;
				v[0]=r[0]+MIN_DISTANCE;
				v[1]=r[2]+MIN_DISTANCE;
			} else if (r[2]==0) {
				v[2]=2;
				v[0]=r[0]+MIN_DISTANCE;
				v[1]=r[1]+MIN_DISTANCE;
			} else {
				std::string s = "GeometryFcc:: vector r has no entry 0: ";
				s = s + PsimagLite::typeToString(r[0]) + " "
					+ PsimagLite::typeToString(r[1]) + " "
					+ PsimagLite::typeToString(r[2]);
				throw std::runtime_error(s.c_str());
			}
		}

		//! output rsum + bsum = input r1 + b1 +r2 +b2
		void add(
				std::vector<RealType> &rsum,
				std::vector<RealType> &bsum,
				const std::vector<size_t>& r1,
				size_t b1,
				const std::vector<size_t>& r2,
				size_t b2) const
		{
			std::vector<RealType> unityvector(3,1);

			std::vector<RealType> b1Plusb2 =
					basisVector_[b1] + basisVector_[b2];
			size_t casetype=computeCase(b1Plusb2);

			switch (casetype) {
			case 0:  // b1+b2 == (a,a,0) or permutations
				rsum =  b1Plusb2 + r1 + r2;
				for (size_t i=0;i<bsum.size();i++) bsum[i] = 0;
				break;
			case 1: // b1 + b2 == (a/2,a/2,a) or permutations
				bsum = unityvector - b1Plusb2;
				rsum =  2.0*b1Plusb2-unityvector + r1 + r2;
				break;
			case 2: // b1 + b2 == (a/2,a/2,0) or permutations
				for (size_t i=0;i<rsum.size();i++) rsum[i] =r1[i]+r2[i];
				bsum = b1Plusb2;
				break;
			}

		}

		size_t computeCase(const std::vector<RealType>& b) const
		{
			size_t flaghalf=0;
			size_t flagzero=0;

			for (size_t i=0;i<b.size();i++) {
				if (b[i]==0.5) flaghalf=1;
				if (b[i]==0) flagzero=1;
			}
			if (flaghalf==0 && flagzero==1) return  0;
			if (flaghalf==1 && flagzero==0) return 1;
			if (flaghalf==1 && flagzero==1) return 2;
			std::string s = "Geometry::computeCase(...) for add: Error\n";
			s += "Exiting at this point " + std::string(__FILE__) + " "
					+ ttos(__LINE__) + "\n";
			throw std::runtime_error(s.c_str());
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

