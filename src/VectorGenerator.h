
#ifndef VECTOR_GENERATOR_H
#define VECTOR_GENERATOR_H
#include <string>

template<typename T>
class VectorGenerator {

public:	
	VectorGenerator(std::vector<T>& data,int localRank) 
		: data_(data), localRank_(localRank_)
	{
	}
	
	void operator()(T& v)
	{
		v = data_[localRank_];
	}

private:
	const std::vector<T>& data_;
	size_t localRank_;

}; // VectorGenerator

#endif