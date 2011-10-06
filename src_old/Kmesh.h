#ifndef K_MESH_H
#define K_MESH_H
#include <vector>

class Kmesh {
public:
	void init(size_t dim,size_t length) 
	{
		 dim_ = dim;
		 length_ = length;
	}
	
	void calcKVector(std::vector<int>& v,size_t kindex) const 
	{
		v.resize(2);
        	v[0] = kindex % length_;
        	v[1] = int(kindex / length_);
	}

	
	size_t length() const { return length_; }
	
private:
	size_t dim_,length_;
};

#endif //K_MESH_H
