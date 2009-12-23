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

extern int spf_entry(int argc,char *argv[],int ,int);

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
	
	myargv = new char*[2];
	myargv[0] = new char[strlen(exename)+1];
	strcpy(myargv[0],argv[0]);
	myargv[1] = new char[strlen(comm1)+1];
	strcpy(myargv[1],argv[1]);
	if (rank==0) {
		cout<<"myargv= "<<myargv[0]<<" "<<myargv[1]<<endl;
	}
	
	spf_entry(2,myargv,rank,GONZA_MPISize);
	
	cout<<"PROCESS "<<rank<<" HAS ENDED\n";
	// All MPI Threads must wait before calling Finalize
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}
	
	
	
