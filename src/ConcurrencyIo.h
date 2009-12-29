#ifndef CONCURRENCY_IO_H
#define CONCURRENCY_IO_H

template <typename ConcurrencyType>
class ConcurrencyIo {
	public:
		template<typename ConcurrencyParam1Type,typename ConcurrencyParam2Type>
		ConcurrencyIo(const ConcurrencyParam1Type& p1,const ConcurrencyParam2Type& p2)
		{
			setPartialCommunicator(p1);
			setPartialCommunicator(p2);
			std::cerr<<"partialComm_.size="<<partialComm_.size()<<"\n";
		}
		
		template<typename T>
		void vectorPrint(const std::vector<T>& v,const std::string& label,std::ofstream& of) 
		{
			std::vector<T> result=v;
			
			std::vector<T> prevResult=v;
			for (size_t i=0;i<partialComm_.size();i++) {
				ConcurrencyType::MPI_Reduce(prevResult,result,ConcurrencyType::MPISUM,0,partialComm_[i]);
				for (size_t j=0;j<result.size();j++) result[j] /= commSize_[i];
				prevResult = result;
			}
			
			for (size_t i=0;i<partialComm_.size();i++) 
				if (rankInComm_[i]!=0) return;
			
			::vectorPrint(result,label,of);
		}
		
		template<typename T>
		void vectorPrint(const T& v,const std::string& label,std::ofstream& of) const
		{
			T result;
			T prevResult = v;
			for (size_t i=0;i<partialComm_.size();i++) {
				ConcurrencyType::MPI_Reduce(&prevResult,&result,1,ConcurrencyType::MPISUM,0,partialComm_[i]);
				result /= commSize_[i];
				prevResult = result;
			}
			for (size_t i=0;i<partialComm_.size();i++) 
				if (rankInComm_[i]!=0) return;
			
			of << label <<" = "<<result << "\n";
		}

	private:
		std::vector<typename ConcurrencyType::MPIComm> partialComm_;
		std::vector<size_t> commSize_;
		std::vector<size_t> rankInComm_;
		
		
		template<typename ConcurrencyParamType>
		void setPartialCommunicator(const ConcurrencyParamType& p)
		{
			if (p.separate()) return;
			std::cerr<<"doing somthhing\n";
			partialComm_.push_back(p.mpiComm());
			commSize_.push_back(p.mpiCommSize());
			rankInComm_.push_back(p.mpiRank());
		}
};

#endif
