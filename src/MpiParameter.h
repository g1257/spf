
template<typename T,typename GeneratorType>
class MpiParameter {

public:
	MpiParameter(T& param,GeneratorType& generator) : generator(param)
	{
	}


}; // MpiParameter

