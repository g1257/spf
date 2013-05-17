
#include "Integrate.h"
#include "ProgramGlobals.h"

template<typename RealType>
struct DenominatorParams
{
	DenominatorParams(const RealType& b,
	                  const RealType& l0,
	                  const RealType& l1)
	    : beta(b),lambda0(l0),lambda1(l1)
	{}

	RealType beta;
	RealType lambda0;
	RealType lambda1;
};

template<typename RealType_>
class DenominatorFunction {

public:

	typedef RealType_ RealType;
	typedef DenominatorParams<RealType> ParametersType;

	DenominatorFunction(ParametersType& params)
	: params_(params)
	{}

	ParametersType& params()
	{
		return params_;
	}

	RealType operator()(RealType x) const
	{
		RealType tmp = 1.0 + exp(-params_.beta * params_.lambda1 * x);
		RealType tmp2 = exp(-params_.beta * params_.lambda0 * x * x);

		return tmp2 * Spf::ProgramGlobals::square(tmp);
	}

private:

	ParametersType& params_;
}; // class ActionFunctor

template<typename RealType_>
class NumeratorFunction {

public:

	typedef RealType_ RealType;
	typedef DenominatorParams<RealType> ParametersType;

	NumeratorFunction(ParametersType& params)
	    : params_(params),den_(params)
	{}

	ParametersType& params()
	{
		return params_;
	}

	RealType operator()(RealType x) const
	{
		return den_(x) * x;
	}

private:

	ParametersType& params_;
	DenominatorFunction<RealType> den_;

}; // class ActionFunctor

int main(int argc,char *argv[])
{
	if (argc<4) {
		std::cout<<"USAGE: "<<argv[0]<<" beta lambda0 lambda1\n";
		return 1;
	}
	DenominatorParams<double> dummy(atof(argv[1]),atof(argv[2]),atof(argv[3]));
	DenominatorFunction<double> function1(dummy);
	double result = Spf::integrate(function1,0.0,2.0);
	std::cout<<result<<"\n";

	NumeratorFunction<double> function2(dummy);
	double result2 = Spf::integrate(function2,0.0,2.0);
	std::cout<<result2<<"\n";
	std::cout<<(result2/result)<<"\n";
}


