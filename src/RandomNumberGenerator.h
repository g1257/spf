#ifndef RNG_HEADER_H
#define RNG_HEADER_H

#include "KISS.h"


class RandomNumberGenerator {
	public:
	void myRandomSeed(unsigned  int seed)
	{
// 		size_t x = 4294967295;
// 		if (rngKiss.getLength()!=x) {
// 			cerr<<"SERIOUS PROBLEM: RANDOM NUMBER GENERATOR NOT RELIABLE ON THIS PLATFORM!!\n";
// 			cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
// 			cerr<<"I will exit now, you should call your customer support 1-800 number.\n";
// 			exit(1);
// 		}
		rngKiss.seed(seed);
		
	}
	
	double myRandom()
	{
		return rngKiss();
	}
	
	private:
	psimag::KISS rngKiss;
	
	
};

#endif
