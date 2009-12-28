#ifndef MPI_IO_H
#define MPI_IO_H

template <typename MpiSystemType>
class MpiIo {
	public:
		template<typename MpiParam1Type,typename MpiParam2Type>
		MpiIo(const MpiParam1Type& p1,const MpiParam2Type& p2) 
		{
			setPartialCommunicator(p1);
			setPartialCommunicator(p2);
		}
		
		template<typename T>
		void vectorPrint(const std::vector<T>& v,const std::string& label,std::ofstream& of)
		{
			std::vector<T> result(v.size());
			std::vector<T> prevResult=v;
					
			for (size_t i=0;i<partialComm_.size();i++) {
				MpiSystemType::MPI_Reduce(&(prevResult[0]),&(result[0]),v.size(),MpiSystemType::MPI_SUM,0,partialComm_[i]);
				for (size_t j=0;j<result.size();j++) result[j] /= commSize_[i];
				prevResult = result;
			}
			
			for (size_t i=0;i<partialComm_.size();i++) 
				if (rankInComm_[i]!=0) return;
			
			vectorPrint(result,label,of);
		}
		
		template<typename T>
		void vectorPrint(const T& v,const std::string& label,std::ofstream& of)
		{
			T result;
			T prevResult = v;
			for (size_t i=0;i<partialComm_.size();i++) {
				MpiSystemType::MPI_Reduce(&prevResult,&result,1,MpiSystemType::MPI_SUM,0,partialComm_[i]);
				result /= commSize_[i];
				prevResult = result;
			}
			for (size_t i=0;i<partialComm_.size();i++) 
				if (rankInComm_[i]!=0) return;
			
			of << label <<" = "<<result << "\n";
		}

	private:
		std::vector<typename MpiSystemType::MPI_Comm> partialComm_;
		std::vector<size_t> commSize_;
		std::vector<size_t> rankInComm_;
		
		template<typename MpiParameterType>
		void setPartialCommunicator(const MpiParameterType& p)
		{
			if (p.separate()) return;
			partialComm_.push_back(p.mpiComm());
			commSize_.push_back(p.mpiCommSize());
			rankInComm_.push_back(p.mpiRank());
		}
};

#endif
