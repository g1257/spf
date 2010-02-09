
/** \ingroup SPF */
/*@{*/

/*! \file RandomNumberGenerator.h
 *
 *  
 *
 */

#ifndef RNG_HEADER_H
#define RNG_HEADER_H

#include "KISS.h"

namespace Spf {
	template<typename FieldType>
	class RandomNumberGenerator {
		public:
		void myRandomSeed(size_t seed)
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
		
		FieldType operator()()
		{
			return rngKiss();
		}
		
		void randomModulus(std::vector<int> &v,size_t conc,int n)
		{
			size_t j;
			v.assign(v.size(),0);
			for (size_t i=0;i<conc;i++) {
				j=0;
				while (v[j]==1) j=size_t(n*this->operator()());
				v[j]=1;
			}
		}
		
		private:
		psimag::KISS rngKiss;
		
		
	};

} // namespace Spf

/*@}*/
#endif
