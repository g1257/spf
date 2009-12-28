#ifndef MPI_SYSTEM_MPI_H
#define MPI_SYSTEM_MPI_H

class MpiSystemMpi {
public:
	typedef ::MPI_Comm MPI_Comm;
	static const int MPI_SUM = ::MPI_SUM;
	static const MPI_Comm MPI_COMM_WORLD = ::MPI_COMM_WORLD;
	
	static int MPI_Reduce(double* operand,double* result,size_t x,size_t op,size_t root ,MPI_Comm comm) 
	{
		MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Reduce(std::complex<double>* operand,std::complex<double>* result,size_t x,size_t op,size_t root ,MPI_Comm comm) 
	{
		MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Reduce(std::vector<double>* operand,std::vector<double>* result,size_t x,size_t op,size_t root ,MPI_Comm comm) 
	{
		MPI_Reduce(operand,result,x,MPI_DOUBLE,x,root,comm);
	}
	
	static int MPI_Comm_Split(MPI_Comm comm1,size_t x,size_t y,MPI_Comm* comm2) 
	{
		MPI_Comm_Split(comm1,x,y,comm2);
	}
};

#endif