
#ifndef CONCURRENCY_PARAM_H
#define CONCURRENCY_PARAM_H

template<typename T,typename GeneratorType,typename ParametersType,typename ConcurrencyType>
class ConcurrencyParameter {

public:
	enum {SEPARATE=0,GATHER=1};
	
	ConcurrencyParameter(T& param,ParametersType& ether,GeneratorType& generator,size_t separateOrTogether,size_t localRank,
		     typename ConcurrencyType::MPIComm& mpiComm,size_t size) 
	: separateOrTogether_(separateOrTogether), localRank_(localRank),mpiCommSize_(size),mpiComm_(mpiComm)
	{
		generator(param);
		if (separateOrTogether_ == SEPARATE) {
			ether.rootname = ether.rootname  + ttos(localRank);
			mpiComm_ = 0;
		}
		
	}
	
	bool separate() const { return (separateOrTogether_==SEPARATE); }
	
	typename ConcurrencyType::MPIComm mpiComm() const { return mpiComm_; }
	
	size_t mpiCommSize() const { return mpiCommSize_; }
	
	size_t mpiRank() const { return localRank_; }

private:
	size_t 	separateOrTogether_,localRank_;
	size_t mpiCommSize_;
	typename ConcurrencyType::MPIComm mpiComm_;
}; // ConcurrencyParameter

#endif

