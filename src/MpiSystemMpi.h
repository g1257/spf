#ifndef MPI_SYSTEM_MPI_H
#define MPI_SYSTEM_MPI_H

#include <mpi.h>

class MpiSystemMpi {
public:
	typedef MPI_Comm MPIComm;
	typedef MPI_Op MPIOp;
	static const MPIOp MPISUM;
	static const MPIComm MPICOMMWORLD;
	
	static int MPI_Reduce(double operand,double& result,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&operand,&result,1,MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Reduce(std::complex<double>& operand,std::complex<double>& result,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&operand,&result,2,MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Reduce(std::vector<double>& operand,std::vector<double>& result,MPI_Op op,size_t root ,MPI_Comm comm) 
	{
		::MPI_Reduce(&(operand[0]),&(result[0]),result.size(),MPI_DOUBLE,op,root,comm);
	}
	
	static int MPI_Reduce(std::vector<std::complex<double> >& operand,std::vector<std::complex<double> >& result,MPI_Op op,size_t root ,MPI_Comm comm)
	{
		size_t x = 2*result.size();
		::MPI_Reduce(&(operand[0]),&(result[0]),x,MPI_DOUBLE,op,root,comm);
	}

	static int MPI_Comm_split(MPIComm comm1,int x,int y,MPIComm* comm2) 
	{
		::MPI_Comm_split(comm1,x,y,comm2);
	}
}; // MpiSystemMpi

const MpiSystemMpi::MPIOp MpiSystemMpi::MPISUM = MPI_SUM;
const MpiSystemMpi::MPIComm MpiSystemMpi::MPICOMMWORLD = MPI_COMM_WORLD;

#endif

