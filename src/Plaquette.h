#ifndef PLAQUETTE_H
#define PLAQUETTE_H
#include "basic.h"
#include <stdexcept>

template<typename DistanceType,typename GeometryType>
class Plaquette {
public:
	//! ugly because geometry needs plaquette
	Plaquette(bool isEnabled) : isEnabled_(isEnabled) { }
	
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
	}
	
	//! Is i in the plaquette given by plaquetteIndex?
	bool isInPlaquette(size_t plaquetteIndex,size_t i) const
	{
		if (!isEnabled_) return false; // should be a warning
		DistanceType rp;
		geometry_->index2Coor(rp,plaquetteIndex,lt_);
		DistanceType ri;
		geometry_->index2Coor(ri,i,lt_);
		return isInPlaquette2(rp,ri);
	}
	
	size_t distance(size_t i,size_t j) const
	{
		if (!isEnabled_) return 0;
		DistanceType ri,rj;
		geometry_->index2Coor(ri,i,lt_);
		geometry_->index2Coor(rj,j,lt_);
		DistanceType dist;
		for (size_t di = 0; di<dim_;di++) {
			int tmp = ri[di]-rj[di];
			//correct for boundary condition;
			if (tmp>=sides_[di]) tmp = latticeLength_[di] - tmp;
			if (tmp<= -sides_[di]) tmp = -latticeLength_[di] - tmp;
			// we add this number so that it is non-negative 
			tmp += sides_[di] - 1;
			// hopefully now tmp is greater than 0
			dist[di] = tmp;
		}
		// The max for $dist[0] is (2*$GlobalLc[0] - 2)
		// and so there are (2*$GlobalLc[0] - 1) of $dist[0]
		return  dist[0] + dist[1]*(2*sides_[0] - 1);
	}
	
	void calcD(size_t j,std::vector<size_t>& d) const
	{
		size_t something  = 2*sides_[0] - 1;
		d[0] = j % something;
		d[1] = size_t(j/something);
		d[0] -= sides_[0] - 1;
		d[1] -= sides_[1] - 1;
		
		
	}
	
private:
	const GeometryType* geometry_;
	bool isEnabled_;
	std::string lt_;
	size_t dim_;
	std::vector<int> latticeLength_;
	std::vector<size_t> sides_;
	
	
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
}; // class Plaquette


#endif
