
#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H
#include <string>
#include <stdexcept>

template<typename FieldType>
class RandomGenerator {

public:	
	RandomGenerator(const std::string& kindOfGenerator,FieldType l0,FieldType l1,int localSize,int localRank) 
		: kindOfGenerator_(kindOfGenerator),l0_(l0), l1_(l1), localRank_(localRank)
	{
		if (kindOfGenerator != "bimodal") throw std::runtime_error("Not supported!!\n");
		std::vector<int> seeds(localSize);
		for (size_t i=0;i<seeds.size();i++) {
			//std::cerr<<"About to do i="<<i<<" total = "<<seeds.size()<<"\n";
			seeds[i] = random();
		}
		//exit(3);
		srand(seeds[localRank]);
	}
	
	void operator()(std::vector<FieldType>& v)
	{
		for (size_t i=0;i<v.size();i++) {
			FieldType x = random()/static_cast<FieldType>(RAND_MAX);
			//std::cerr<<"x="<<x<<" "<<l0_<<" "<<l1_<<"\n";
			if (x<0.5) v[i] = l0_ + l1_;
			else v[i] = l0_ - l1_;
		}
	}

private:
	std::string kindOfGenerator_;
	FieldType l0_,l1_,localRank_;

}; // RandomGenerator

#endif

