#ifndef CONCURRENCY_SERIAL_H
#define CONCURRENCY_SERIAL_H

class ConcurrencySerial {
	typedef std::complex<double> ComplexType;
	
public:
	typedef  int MPIComm;
	static const int MPISUM = 0;
	static const int MPICOMMWORLD = 0;
	
	static int MPI_Reduce(double operand,double& result,int op,size_t root,MPIComm comm) 
	{
		//MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Reduce(std::complex<double> operand,std::complex<double>& result,int op,size_t root ,MPIComm comm) 
	{
		//MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Reduce(std::vector<double>& operand,std::vector<double>& result,int op,size_t root ,MPIComm comm) 
	{
		//MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Reduce(std::vector<ComplexType>& operand,std::vector<ComplexType>& result,int op,size_t root ,MPIComm comm) 
	{
		//MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Comm_split(MPIComm comm1,size_t x,size_t y,MPIComm* comm2) 
	{
		//MPI_Comm_Split(comm1,x,y,comm2);
	}
	
};

#endif
