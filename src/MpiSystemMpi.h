#ifndef MPI_SYSTEM_MPI_H
#define MPI_SYSTEM_MPI_H

#include <mpi.h>

class MpiSystemMpi {
public:
	typedef MPI_Comm MPIComm;
	typedef MPI_Op MPIOp;
	static const MPIOp MPI_SUM;
	static const MPIComm MPI_COMM_WORLD;
	
	static int MPI_Reduce(double operand,double& result,int x,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&operand,&result,x,MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Reduce(std::complex<double>& operand,std::complex<double>& result,size_t x,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&operand,&result,x,MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Reduce(std::vector<double>& operand,std::vector<double>& result,size_t x,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&(operand[0]),&(result[0]),x,MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Comm_split(MPIComm comm1,int x,int y,MPIComm* comm2) 
	{
		::MPI_Comm_split(comm1,x,y,comm2);
	}
}; // MpiSystemMpi

const MpiSystemMpi::MPIOp MpiSystemMpi::MPISUM = MPI_SUM;
//const MPIOp MpiSystemMpi::MPI_SUM = MPI_SUM;
//MpiSystemMpi::MPIComm MpiSystemMpi::MPI_COMM_WORLD = MPI_COMM_WORLD;

#endif

