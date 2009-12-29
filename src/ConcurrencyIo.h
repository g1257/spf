#ifndef CONCURRENCY_IO_H
#define CONCURRENCY_IO_H
#include <mpi.h>

template <typename ConcurrencyType_>
class ConcurrencyIo {
	public:
		typedef ConcurrencyType_ ConcurrencyType;
	
		template<typename ConcurrencyParam1Type,typename ConcurrencyParam2Type>
		ConcurrencyIo(const ConcurrencyParam1Type& p1,const ConcurrencyParam2Type& p2)
		{
			setPartialCommunicator(p1);
			setPartialCommunicator(p2);
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			int color = int(rank / 3);
			MPI_Comm_split(MPI_COMM_WORLD,color,rank,&testComm_);
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
				ConcurrencyType::MPI_Reduce(prevResult,result,ConcurrencyType::MPISUM,0,partialComm_[i]);
				result /= commSize_[i];
				prevResult = result;
			}
			for (size_t i=0;i<partialComm_.size();i++) 
				if (rankInComm_[i]!=0) return;
			
			of << label <<" = "<<result << "\n";
		}

		bool canWrite() const
		{
			for (size_t i=0;i<partialComm_.size();i++)
                                if (rankInComm_[i]!=0) return false;
			
			return true;
		}

		template<typename SomeType>
		bool average(SomeType& value)
		{
			SomeType result = value;
                        //SomeType prevResult = value;
			size_t i = 0;
                        //for (size_t i=0;i<partialComm_.size();i++) {
				ConcurrencyType::MPI_Reduce(value,result,ConcurrencyType::MPISUM,0,testComm_);
				//result /= commSize_[i];
				//prevResult = result;
			//}
			value = result;
			
			for (size_t i=0;i<partialComm_.size();i++)
                                if (rankInComm_[i]!=0) return false;
			
			return true;

		}

		bool average(double& value)
		{
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			std::cerr<<"We're here "<<rank<<" "<<value<<"\n";
			double result = value;
			MPI_Reduce(&value,&result,1,MPI_DOUBLE,MPI_SUM,0,testComm_);
			std::cerr<<"Average "<<rank<<" "<<result<<"\n";
			value = result;
		}

	private:
		std::vector<typename ConcurrencyType::MPIComm> partialComm_;
		std::vector<size_t> commSize_;
		std::vector<size_t> rankInComm_;
		MPI_Comm testComm_;
		
		template<typename ConcurrencyParamType>
		void setPartialCommunicator(const ConcurrencyParamType& p)
		{
			if (p.separate()) return;
			//std::cerr<<"doing somthhing\n";
			partialComm_.push_back(p.mpiComm());
			commSize_.push_back(p.mpiCommSize());
			rankInComm_.push_back(p.mpiRank());
		}
};

#endif
