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





#ifndef IO_H_
#define IO_H_

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <functional>
#include <cstring>
#include "basic.h"
#include "parameters.h"
#include "dynvars.h"
#include "aux.h"
#include "geometry.h"
#include "argument.h"

#ifdef __LIBCATAMOUNT__ 
#include <catamount/dclock.h>
#endif

/*! \brief A class for the handling of input/output 

 * This is a singleton, that is a class that can
 * be instantiated only once. Always pass it as a reference. */

class Io {
	public:
			Io();
			
			~Io();
	
			void initOutput(Parameters &ether);
			
			void printMiscInfo()
			{
				getHostInfo();
				getCompilerInfo();
				for (int i=0;i<nFiles;i++) {
					if (isInVector(i,excludedfiles)) continue;
					printMiscInfo(file[i]);
				}
			}
			
			void printMiscInfo(ostream &s);
			
			
			void StartTimer()
			{
				time_t *tp = new time_t;
				tref = time(tp);
				delete tp;
			}
			
			void getCompilerInfo() 
			{
				strcpy(compInfo,__DATE__);
				strcat(compInfo," ");
				strcat(compInfo,__TIME__);
				strcat(compInfo," ");
			}
			
			void writeFinalStuff();
		
			void currentTime();
		
			ostream &fileError() const
			{
#ifndef _AIX
				// NullStream ns;
				// if (rank>0) return ns;
#endif
				return std::cerr;
			}
			
			ostream &fileOut() const
			{
#ifndef _AIX
				// NullStream ns;
				// if (rank>0) return ns;
#endif
				return std::cout;
			}
			
			void historyPrint(string const &s,int option=0)
			{
				if (!doesPrintHistory) return;
				file[0]<<s;
				if (option==0) file[0]<<endl;
			}
			
			
			void printAverages(Parameters &ether,Aux &aux);
			
			void printSnapshot(DynVars const &dynVars,Parameters const &ether);
			void printSnapshot(DynVars const &dynVars,Parameters const &ether,int option);
			void printSnapshot(DynVars const &dynVars,Parameters const &ether,std::ofstream  &f);
			void printLcd(Parameters const &ether,Aux &aux);
			void getHostInfo();
			int input(char const *filename,Geometry &geometry,DynVars &dynVars,Parameters &ether,Aux &aux);
			
	private:
			int setBoundaryConditions(std::string const &s,Parameters &ether);
			std::ofstream *file;
			int nFiles;
			char filename[256],hostname[256],compInfo[256],cTime[256];
			time_t tref;
			static bool isInstatiated,isInit;
			std::string *extensions;
			int rank,rankTpem;
			bool doesPrintHistory;
			std::vector<int> excludedfiles;
};

void MPI_VecReduce(Histogram &v1,Histogram &v2,int n,int type,int op,int
rank,int comm);

#endif
