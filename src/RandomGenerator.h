
#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H
#include <string>

class RandomGenerator {

public:	
	RandomGenerator(const std::string& kindOfGenerator,size_t l0,size_t l1,int localSize,int localRank) 
		: kindOfGenerator_(kindOfGenerator),l0_(l0), l1_(l1), localRank_(localRank_)
	{
		if (kindOfGenerator != "bimodal") throw std::runtime_error("Not supported!!\n");
		std::vector<int> seeds(localSize);
		for (size_t i=0;seeds.size();i++) seeds[i] = random();
		srand(seeds[localRank]);
	}
	
	void operator()(std::vector<double>& v)
	{
		for (size_t i=0;i<v.size();i++) {
			if (random()<0.5) v[i] = l0_ + l1_;
			else v[i] = l0_ - l1_;
		}
	}

private:
	std::string kindOfGenerator_;
	size_t l0_,l1_,localRank_;

}; // RandomGenerator

#endif

