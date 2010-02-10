
/** \ingroup SPF */
/*@{*/

/*! \file GreenFunction.h
 *
 *  FIXME
 *
 */
#ifndef GREENFUNCTION_H
#define GREENFUNCTION_H
#include "Utils.h"

namespace Spf {
	template<typename AlgorithmType>
	class GreenFunction {
		typedef typename AlgorithmType::FieldType FieldType;
		typedef typename AlgorithmType::ComplexType ComplexType;
	public:	
		
		GreenFunction(AlgorithmType& algorithm) : algorithm_(algorithm)
		{
		}
		
		ComplexType operator()(size_t lambda1,size_t lambda2)
		{
			return algorithm_.greenFunction(lambda1,lambda2);
		}
		
	private:			
		AlgorithmType& algorithm_;
		
	}; // GreenFunction
	
	template<typename AlgorithmType>
	std::ostream& operator<<(std::ostream& os,GreenFunction<AlgorithmType>& algorithm)
	{
		throw std::runtime_error("unimplemented operator<< for GreenFunction\n");
		return os;
	}
} // namespace Spf

/*@}*/
#endif
