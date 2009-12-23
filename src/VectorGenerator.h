
#ifndef VECTOR_GENERATOR_H
#define VECTOR_GENERATOR_H
#include <string>

template<typename T>
class VectorGenerator {

public:	
	VectorGenerator(std::vector<T>& data,size_t localRank) 
		: data_(data), localRank_(localRank)
	{
		//std::cerr<<"localRank="<<localRank_<<"\n";
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