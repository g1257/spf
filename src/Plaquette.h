#ifndef PLAQUETTE_H
#define PLAQUETTE_H
#include "basic.h"
#include "Matrix.h"
#include <stdexcept>

template<typename DistanceType,typename GeometryType>
class Plaquette {
public:
	//! ugly because geometry needs plaquette
	Plaquette(bool isEnabled) : isEnabled_(isEnabled)
	{
	}
	
	void init(const GeometryType* geometry,size_t side)
	{
		if (side==0) isEnabled_ = false;
		if (!isEnabled_) return;
		geometry_ = geometry;
		lt_ = geometry->latticeType;
		dim_ = geometry->dim();
		latticeLength_ = geometry->Length;
		sides_.clear();
		for (size_t i=0;i<dim_;i++) sides_.push_back(side);
		
		if (lt_!="square") throw std::runtime_error("Plaquette unimplemented\n");
		size_t n = geometry->volume();
		isInPlaquetteCached_.resize(n,n);	
                for (size_t plaquetteIndex=0;plaquetteIndex<n;plaquetteIndex++)
                        for (size_t i=0;i<n;i++)
                                isInPlaquetteCached_(plaquetteIndex,i) = isInPlaquette_(plaquetteIndex,i);

	}
	
	//! Is i in the plaquette given by plaquetteIndex?
	bool isInPlaquette(size_t plaquetteIndex,size_t i) const
	{
		if (!isEnabled_) return false; // should be a warning
		return isInPlaquetteCached_(plaquetteIndex,i);
	}
	
	size_t distance(size_t i,size_t j) const
	{
		if (!isEnabled_) return 0;
		DistanceType ri,rj;
		geometry_->index2Coor(ri,i,lt_);
		geometry_->index2Coor(rj,j,lt_);
		std::vector<int> dist(dim_);
		for (size_t di = 0; di<dim_;di++) {
			int tmp = ri[di]-rj[di];
			//correct for boundary condition;
			if (tmp>=sides_[di]) tmp = latticeLength_[di] - tmp;
			int minusSide = sides_[di];
			minusSide *= (-1);
			if (tmp<= minusSide) tmp = -latticeLength_[di] - tmp;
			dist[di] = tmp;
			if (dist[di] >= sides_[di] || dist[di]<=minusSide) throw std::runtime_error("distance\n");
			
		}
		// The max for $dist[0] is (2*$GlobalLc[0] - 2)
		// and so there are (2*$GlobalLc[0] - 1) of $dist[0]
		return packDistance(dist); 
	}
	
	void unpackDistance(std::vector<int>& d,size_t distanceIndex) const
	{
		size_t something  = 2*sides_[0] - 1;
		d.resize(2);
		d[0] = distanceIndex % something;
		d[1] = int(distanceIndex/something);
		d[0] -= (sides_[0] - 1);
		d[1] -= (sides_[1] - 1);
		
		
	}
	
	size_t packDistance(std::vector<int>& d) const
	{
		// we add this number so that it is non-negative 
		d[0] += (sides_[0] - 1);
		d[1] += (sides_[1] - 1);
		
		size_t something  = 2*sides_[0] - 1;
		
		int x =  d[0] + d[1]*something;
		if (x<0) throw std::runtime_error("packDistance\n");
		return x;
	}
	
private:
	const GeometryType* geometry_;
	bool isEnabled_;
	std::string lt_;
	size_t dim_;
	std::vector<int> latticeLength_;
	std::vector<int> sides_;
	psimag::Matrix<bool> isInPlaquetteCached_;	
	
	
	bool isInPlaquette2(const DistanceType& rp,const DistanceType& ri) const
	{
		for (size_t di = 0; di<dim_;di++) 
			if (!isInPlaquette3(rp,ri,di)) return false ;
		
		return true;
	}

	bool isInPlaquette3(const DistanceType& rp,const DistanceType& ri,size_t whatDimension) const
	{
		size_t tmp = vecDistOneDim(ri,rp,whatDimension); // wraps around
		if (tmp>=sides_[whatDimension]) return false ;
		return true;
	}
	
	size_t vecDistOneDim(const DistanceType& ri,const DistanceType& rj,size_t whatDimension) const
	{
		int tmp = ri[whatDimension]-rj[whatDimension];
		if (tmp<0) tmp += latticeLength_[whatDimension];
		return tmp;
	}
	
	//! Is i in the plaquette given by plaquetteIndex?
        bool isInPlaquette_(size_t plaquetteIndex,size_t i) const
        {
                if (!isEnabled_) return false; // should be a warning
                DistanceType rp;
                geometry_->index2Coor(rp,plaquetteIndex,lt_);
                DistanceType ri;
                geometry_->index2Coor(ri,i,lt_);
                return isInPlaquette2(rp,ri);
        }

}; // class Plaquette


#endif
