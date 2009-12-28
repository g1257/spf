
#ifndef MPI_PARAMETER_H
#define MPI_PARAMETER_H

template<typename T,typename GeneratorType,typename ParametersType,typename MpiSystemType>
class MpiParameter {

public:
	enum {SEPARATE,TOGETHER};
	
	MpiParameter(T& param,ParametersType& ether,GeneratorType& generator,size_t separateOrTogether,size_t localRank,
		     typename MpiSystemType::MPI_Comm& mpiComm,size_t size) 
	: separateOrTogether_(separateOrTogether), localRank_(localRank),mpiCommSize_(size),mpiComm_(mpiComm)
	{
		generator(param);
		if (separateOrTogether_ == SEPARATE) {
			ether.rootname = ether.rootname  + ttos(localRank);
			mpiComm_ = 0;
		}
		
	}
	
	bool separate() const { return (separateOrTogether_==SEPARATE); }
	
	typename MpiSystemType::MPI_Comm mpiComm() const { return mpiComm_; }
	
	size_t mpiCommSize() const { return mpiCommSize_; }
	
	size_t mpiRank() const { return localRank_; }

private:
	size_t 	separateOrTogether_,localRank_;
	size_t mpiCommSize_;
	typename MpiSystemType::MPI_Comm mpiComm_;
}; // MpiParameter

#endif

