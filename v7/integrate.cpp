
#include "Integrate.h"

class ActionParams
{};

template<typename RealType_>
class ActionFunctor {

public:

	typedef RealType_ RealType;
	typedef ActionParams ParametersType;

	ActionFunctor(ActionParams& params)
	: params_(params)
	{}

	ActionParams& params()
	{
		return params_;
	}

	RealType operator()(RealType x) const
	{
		return x;
	}

private:

	ActionParams& params_;
}; // class ActionFunctor

int main(int argc,char *argv[])
{
	ActionParams dummy;
	ActionFunctor<double> function1(dummy);
	double result = Spf::integrate(function1,0.0,2.0);
	std::cout<<result<<"\n";
}


