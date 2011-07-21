
/** \ingroup SPF */
/*@{*/

/*! \file RootFindingBisection.h
 *
 *  RootFindingBisection such as chemical potential
 *
 */
#ifndef ROOT_FIND_BISECT_H
#define ROOT_FIND_BISECT_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite

namespace Spf {
	template<typename FunctionType>
	class RootFindingBisection {
		typedef typename FunctionType::FieldType FieldType;
		
	public:
		RootFindingBisection(const FunctionType& function,
		                  size_t maxIter=100,
		                  FieldType tolerance=1.0e-3)
		: function_(function),maxIter_(maxIter),tolerance_(tolerance),
		  a_(-100),b_(100)
		{
		}
		
		void operator()(FieldType& mu)
		{
			// INPUT: Function f, endpoint values a, b, tolerance TOL, 
			//                maximum iterations NMAX
			// CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or 
			//                           f(a) > 0 and f(b) < 0
			// OUTPUT: value which differs from a root of f(x)=0 by less than TOL
		
			size_t n = 0;
			while(n<maxIter_) {
				FieldType c = (a_ + b_)/2; // new midpoint
				
				if (function_(c) == 0 || (b_ - a_)/2 < tolerance_) { //solution found
					mu = c;
					return;
				}
				n++;
				if (sign(function_(c)) == sign(function_(a_)))
					a_ = c;
				else
					b_ = c;
			}
			throw std::runtime_error("RootFindBisection failed\n");
		}
		
	private:
		
		int sign(const FieldType& x) const
		{
			return (x>=0) ?  1  : -1;
		}
		
		const FunctionType& function_;
		size_t maxIter_;
		FieldType tolerance_;
		FieldType a_,b_;
	}; // RootFindingBisection
} // namespace Spf

/*@}*/
#endif// ROOT_FIND_BISECT_H

