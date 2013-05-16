
/** \ingroup SPF */
/*@{*/

/*! \file AlgorithmFactory.h
 *
 *  Monte Carlo for SPF
 *
 */
#ifndef ALGORITHM_FACTORY_H
#define ALGORITHM_FACTORY_H
#include "ProgressIndicator.h"
#include <cassert>

namespace Spf {
	template<typename Type1,typename Type2>
	class AlgorithmFactory {

	public:

		AlgorithmFactory(Type1* t1,Type2* t2) : t1_(t1), t2_(t2)
		{
			assert((t1 && !t2) || (t2 && !t1));
		} 

		void init()
		{
			(t1_) ? t1_->algorithm().init() : t2_->algorithm().init();
		}

		template<typename SomeRngType,typename SomeOperationsType>
		bool isAccepted(size_t i,SomeRngType& rng,SomeOperationsType& ops,int n)
		{
			return (t1_) ? t1_->algorithm().isAccepted(i,rng,ops,n)
			             : t2_->algorithm().isAccepted(i,rng,ops,n);
		}

		template<typename SomeOperationsType>
		void accept(size_t i,SomeOperationsType& ops)
		{
			return (t1_) ? t1_->algorithm().accept(i,ops)
			             : t2_->algorithm().accept(i,ops);
		}

	private:

		Type1* t1_;
		Type2* t2_;
	}; // class AlgorithmFactory
} // namespace Spf

/*@}*/
#endif // ALGORITHM_FACTORY_H

