
#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H
#include <string>
#include <stdexcept>

template<typename FieldType>
class RandomGenerator {
	static const size_t MAX_SEEDS = 1000;
public:	
	RandomGenerator(const std::string& kindOfGenerator)  
		: kindOfGenerator_(kindOfGenerator),seeds_(MAX_SEEDS)
	{
		if (kindOfGenerator != "bimodal") throw std::runtime_error("Not supported!!\n");
		for (size_t i=0;i<seeds_.size();i++) {
			//std::cerr<<"About to do i="<<i<<" total = "<<seeds.size()<<"\n";
			seeds_[i] = random();
		}
	}

	void setBimodal(FieldType l0,FieldType l1) 
	{
		l0_=l0;
		l1_=l1;
	}
	
	void seed(size_t number) 
	{
		srand(seeds_[number]);
		std::cerr<<"I'm "<<number<<"with seed="<<seeds_[number]<<"\n";
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
	std::vector<int> seeds_;
	FieldType l0_,l1_;
	

}; // RandomGenerator

#endif

