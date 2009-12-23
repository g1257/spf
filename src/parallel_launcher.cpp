/*******************************************************************
 **         THIS FILE IS PART OF THE SPF (SPIN-PHONON-FERMION)    ** 
 **                          COMPUTER PROGRAM                     **
 **                                                               **
 **			 SPF VERSION 6.4                          **
 **                    Dated: March 29, 2006                      **
 **                                                               **
 **   
                                                             **
FOR INTERNAL USE ONLY. 
YOU MAY NOT DISTRIBUTE OR RE-DISTRIBUTE THIS SOFTWARE IN
ANY FORM.
      
DISCLAIMER

THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS "AS IS" AND ANY 
EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY 
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Note: Please acknowledge SPF if you publish results obtained
with it. You can use the following acknowledgement line:
"This research used the SPF 
computer code (http://mri-fre.ornl.gov/spf)." 
 **                                                               **
 **       For more info see: http://mri-fre.ornl.gov/spf          **
 ******************************************************************/




// #
// # MPI wrapper to launch many serial jobs  as one paralell job
// # by G. Alvarez 

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cstdlib>

extern int spf_entry(int argc,char *argv[]);

#define MAX_WORD 256
using namespace std;

extern void num2digits(char *temp,int &n,int x);

// MAIN PROGRAM BEGINS HERE
int main(int argc,char **argv) 
{
	register int i;
	int myOffset,number;
	int rank,GONZA_MPISize,serial_tasks;
	char exename[MAX_WORD], tempc[MAX_WORD],let[3],comm1[MAX_WORD];
	char **myargv;
	// Since the slave processes get wild from the beginning it
	// doesnt matter where you call Init(...)
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&GONZA_MPISize);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	ether.mpiSize = GONZA_MPISize;
	
	cout<<"Here  "<<rank<<endl;	
	// MPI is so stupid that does not do some trivial things
	// like passing to each spawn process the command line arguments
	// or the current working directory. Those things must be done
	// by hand. Broadcast is used as the communication function.
	// Rank 0 process should get everything correct and is themaster here.
	// First, let's deal with the current working dir
	getcwd(tempc,MAX_WORD);
	MPI_Bcast(tempc,strlen(tempc)+1,MPI_CHAR,0,MPI_COMM_WORLD);
	if (rank!=0) chdir(tempc);
	
	// Ok, now to the command line arguments
	for (i=0;i<argc;i++) {
		MPI_Bcast(argv[i],strlen(argv[i])+1,MPI_CHAR,0,MPI_COMM_WORLD);
	}
	
	if (argc==1) myOffset=0;
	else myOffset=atoi(argv[1]);
        if (argc==3) strcpy(exename,argv[2]);
	else strcpy(exename,"orb1new");
 
	number=rank+myOffset;	
	int nn;
	num2digits(let,nn,number);
	
	strcpy(comm1,"test");
	strcat(comm1,let);
	strcat(comm1,".inp");
	myargv = new char*[2];
	myargv[0] = new char[strlen(exename)+1];
	strcpy(myargv[0],exename);
	myargv[1] = new char[strlen(comm1)+1];
	strcpy(myargv[1],comm1);
	if (rank==0) {
		cout<<"myargv= "<<myargv[0]<<" "<<myargv[1]<<endl;
	}
	
	spf_entry(2,myargv);
	
	cout<<"PROCESS "<<rank<<" HAS ENDED\n";
	// All MPI Threads must wait before calling Finalize
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}
	
	
	
