
template<typename T,typename GeneratorType,typename ParametersType>
class MpiParameter {

public:
	enum {SEPARATE,TOGETHER};
	
	MpiParameter(T& param,ParametersType& ether,GeneratorType& generator,size_t separateOrTogether,size_t localRank) 
	{
		generator(param);
		if (separateOrTogether == SEPARATE) {
			ether.rootname = ether.rootname  + ttos(localRank);
		}
	}
	
	


}; // MpiParameter

