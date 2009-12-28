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
template<typename MpiIoType>
class Io {
	public:
			
			Io(MpiIoType* mpiIo);
			
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
			MpiIoType* mpiIo_;
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


template<typename MpiIoType>
Io<MpiIoType>::Io(MpiIoType* mpiIo) : mpiIo_(mpiIo)
{
	
	if (isInstatiated) {
		cerr<<"FATAL: Only one Io object per program please!\n";
		cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<endl;
		exit(1);
	}
	nFiles = 11;
	/*! <h2> OUTPUT FILES </h2><br>
	 * STDOUT, STDERR and other files whose names start with the value of the parameter 
	 * OUTPUT_ROOTNAME are written.
	 * To each file a signature, data and footer is written in that order.<br>
	 * Sometimes what is printed depends on the "diagonalization" algorithm used, DIAG or TPEM
	 * as indicated.<br>
	 * <b>Description of the Output Files</b>
	 * -  \b  STDOUT: What is written: Progress
	 * - \b STDERR: What is written: Warnings, errors and debugging information
	 * - \b .dat:  What is written: History.
	 * Use average.pl to get the Averages. For example:
	 * perl average.pl OUTPUT_ROOTNAME.dat "Mag="
	 * Note that if the observable name has parenthesis, brackets or any special characters it 
	  * must be protected. For example to get the structure factor at q=0, do:
	 * perl average.pl OUTPUT_ROOTNAME.dat "S\(0\)="
	 * - \b .now:  What is written:   N(w) for every  independent band if DIAG method was used 
	 * - \b .mom: What is written: The average moments of the TPEM expansion if TPEM was used. 
	 *         Helpful to calculate the DOS using the script mom2nw.pl
	 * - \b .sav:  What is written: The First and Last configuration (dynamical variables)
	 * - \b .lcd:  What is written: Local Charge Density. FIXME: Not Implemented.
	 * - \b .eig:  What is written: Mean value of the eigenvalues.
	 * - \b .cor:  What is written:   The correlation lattice for classical spins.
	        Check source code for normalization. Remember that not all sites need have spins. 
	 * - \b .arw_: A(r,w) under DIAG, off-diagonal moments under TPEM
	 * - \b .ldos: Local Density of States (DIAG Only)
	 * - \b .opt: Optical conductivity (DIAG) or off-diagonal moments to calculate it (TPEM)
	 * - \b .cco: Charge correlations (DIAG only) and only if options do not have "fileminizize" and do have "chargecorrelation"
	 * - \b .oco: Orbital correlations (DIAG only), use with option "orbitalcorrelation"
	 */
	extensions = new std::string[nFiles];
	extensions[0]=".dat"; extensions[1]=".now",extensions[2]=".sav";
	extensions[3]=".arw_";extensions[4]=".ldos"; extensions[5]=".opt";
	extensions[6]=".oco"; extensions[7]=".lcd"; extensions[8]=".eig";
	extensions[9]=".cco"; extensions[10]=".cor";
			
	isInstatiated=true;
	file = new std::ofstream[nFiles];
	StartTimer();
	currentTime();
}

template<typename MpiIoType>
Io<MpiIoType>::~Io()
{
	
	if (isInit) {
		writeFinalStuff();
		for (int i=0;i<nFiles;i++) {
			if (isInVector(i,excludedfiles)) continue;
			if (i>0 && rank>0) break;
			file[i].close();
		}
	}
	if (isInstatiated) {
		delete [] file;
		delete [] extensions;
	}
	isInit=isInstatiated=false;
				
}

template<typename MpiIoType>
void Io<MpiIoType>::initOutput(Parameters &ether)
{
	// Prepare files
	int i;
	doesPrintHistory=true;
	//ether.setOption("fileminimize");
	if (ether.isSet("fileminimize")) {
		excludedfiles.clear();
		excludedfiles.push_back(1);
		excludedfiles.push_back(2);
		excludedfiles.push_back(3);
		excludedfiles.push_back(4);
		excludedfiles.push_back(5);
		excludedfiles.push_back(6);
		excludedfiles.push_back(7);
	}
	
	rank = ether.mpiRank;
	rankTpem=ether.mpiTpemRank;
	if (ether.mpiNop2==1 && rank>0) doesPrintHistory=false;
	if (ether.mpiNop2>1 && ether.mpiTpemRank>0)  doesPrintHistory=false;
	
	if (ether.mpiTpemRank==0) {			
	getHostInfo();
	getCompilerInfo();
				
	std::string filename;
	for (i=0;i<nFiles;i++) {
		//if (i>0 && rank>0) break;
		if (isInVector(i,excludedfiles)) continue;
		filename = std::string(ether.rootname) + extensions[i];
		file[i].open(filename.c_str());
		ether.print(file[i]);
		printMiscInfo(file[i]);
	}
	if (ether.isSet("verbose")) std::cerr<<"Io::init: success, rank="<<rank<<"\n";
	isInit=true;
	}
			
}

template<typename MpiIoType>
void Io<MpiIoType>::printMiscInfo(ostream &s)
{
	int pid=getpid();
	/*! <b>Signature of the output files</b><br>
	 * The input parameters followed by:
	 * - \b #InitialTime: The Time/Date the program started.
	 * - \b #PID: PID for the program (in parallel programs, pid of the master process)
	 * - \b #HOSTNAME Hostname the program runs on. (output of uname -a)
	 * (in parallel programs, hostname where the  master process runs)
	 */

	s<<"#InitialTime="<<cTime<<endl;
	s<<"#PID="<<pid<<endl;
	s<<"#HOSTNAME="<<hostname<<endl;
	s<<"#COMPILED="<<compInfo<<endl;		
}
			
template<typename MpiIoType>
void Io<MpiIoType>::writeFinalStuff()
{
	
        if (rankTpem!=0) return;
        double tuser,tsystem,treal;
        struct tms *buffer = new tms;
        time_t *tp = new time_t;
        long clk_tck= sysconf(_SC_CLK_TCK);
        clock_t t=0;
        time_t stw_t2 = time(tp);
        treal=difftime(stw_t2,tref);

#ifndef __LIBCATAMOUNT__ 
        t=times(buffer);
        tuser=(double)buffer->tms_utime/clk_tck;
        tsystem=(double)buffer->tms_stime/clk_tck;
#else
	t=dclock();
        tuser=tsystem=-1;
#endif

        currentTime();
        /*! <b> Footer printed to all output files</b>:
         * User, System and Real times for the program execution in seconds. */
        for (int i=0;i<nFiles;i++) {
		if (isInVector(i,excludedfiles)) continue;
                file[i]<<"# FinalTime "<<cTime<<endl;
                if (tuser<0 || tsystem<0) {
file[i]<<"# This platform does not support user/system time\n";
file[i]<<"# Read the real timer instead.\n";
                } else {
                        file[i]<<"# User   "<<tuser<<endl;
                        file[i]<<"# System "<<tsystem<<endl;
                }
		file[i]<<"# Real   "<<treal<<endl;
                file[i]<<"# EOF"<<endl;
        }
        delete buffer;
	delete  tp;
}


template<typename MpiIoType>
void Io<MpiIoType>::currentTime()
{
	time_t *tp = new time_t;
	time(tp);
	if (!tp) strcpy(cTime,"ERROR\n");
	else strcpy(cTime,asctime(localtime(tp)));
	int l=strlen(cTime);
	if (l>0) cTime[l-1]='\0';
	delete tp;
}

template<typename MpiIoType>
void Io<MpiIoType>::printLcd(Parameters const &ether,Aux &aux)
{
	vectorPrint(aux.lcd,"lcd",file[7]);
}

template<typename MpiIoType>
void Io<MpiIoType>::printAverages(Parameters &ether,Aux &aux)
{
	int i,j,bandIndex=0,temp;
	
#ifdef USE_MPI
	if (ether.mpiNop2>1) {
		// MPI_VecReduce(aux.lcd); 
		// MPI_VecReduce(aux.clasCor);
		// MPI_VecReduce(aux.avMoments);
		//MPI_VecReduce(aux.Nw[0],aux.Nw[0],aux.Nw[0].size(),-1,-1,0,ether.mpiCommAdj);
	}
	if (ether.mpiTpemRank>0) return;
	
#endif	
	if (!ether.isSet("savelcd")) {	
		vectorDivide(aux.lcd,ether.iterEffective*ether.mpiNop2);
		vectorPrint(aux.lcd,"lcd",file[7]);
		//file[0]<<"#SUMLCD="<<vectorSum(aux.lcd);
	}
	vectorDivide(aux.eigM,ether.iterEffective*ether.mpiNop2);
	vectorPrint(aux.eigM,"eigenvals",file[8]);
	vectorDivide(aux.clasCor,ether.iterEffective*ether.mpiNop2);
	vectorPrint(aux.clasCor,"clasCor",file[10]);
	if (aux.orbitalAngles.size()>0) {
		vectorDivide(aux.orbitalAngles,ether.iterEffective*ether.mpiNop2);
		mpiIo_->vectorPrint(aux.orbitalAngles,"OrbitalAngles",file[0]);
	}
	vectorDivide(aux.avMoments,ether.iterEffective*ether.mpiNop2);
	mpiIo_->vectorPrint(aux.avMoments,"moments",file[0]);
	
	if (ether.bcsDelta0>0) 	{	
		vectorDivide(aux.bcsCorxx,ether.iterEffective*ether.mpiNop2*ether.linSize);
		mpiIo_->vectorPrint(aux.bcsCorxx,"bcscorxx",file[0]);
		
		vectorDivide(aux.bcsCorxy,ether.iterEffective*ether.mpiNop2*ether.linSize);
		mpiIo_->vectorPrint(aux.bcsCorxy,"bcscorxy",file[0]);
	}
	
	for (bandIndex=0;bandIndex<ether.numberOfOrbitals;bandIndex++) {
		file[0]<<"# Band "<<bandIndex<<"\n";
		aux.Nw[bandIndex].divide(ether.iterEffective*ether.mpiNop2,1);
		aux.Nw[bandIndex].print(file[0]);
	}
	
	if (ether.isSet("optical")) {
		if (ether.tpem) {
			vectorDivide(aux.opticalMoments,ether.iterEffective);
			mpiIo_->vectorPrint(aux.opticalMoments,"optMoments",file[5]);
		} else {
			aux.Sigma.divide(ether.iterEffective,1);
			aux.Sigma.print(file[5]);
			cerr<<"printed sigma\n";
		}
	} // optical
	
	
	if (ether.isSet("akw")) {
		if (ether.tpem) {	 
			vectorDivide(aux.offdMoments,ether.iterEffective);
			vectorPrint(aux.offdMoments,"moments",file[3]);
		} else {
			temp=ether.linSize;
			if (ether.typeofmodel=="MODEL_KONDO_DMS_FCC") temp=6*ether.linSize;
			for (i=0;i<temp;i++) {
				aux.Arw[i].divide(ether.iterEffective,1);
				aux.ArwC[i].divide(ether.iterEffective,1);
			}
			for (i=0;i<aux.Arw[0].size();i++)  {
				file[3]<<aux.Arw[0].coorX(i)<<" ";
				for (j=0;j<temp;j++) 		
					file[3]<<aux.Arw[j].coorY(i)<<" ";
				file[3]<<endl;
			}
			file[3]<<"% Imaginary part below\n";
			for (i=0;i<aux.ArwC[0].size();i++)  {
				file[3]<<aux.ArwC[0].coorX(i)<<" ";
				for (j=0;j<temp;j++) 		
					file[3]<<aux.ArwC[j].coorY(i)<<" ";
				file[3]<<endl;
			}
		}
	} // akw
	
	if (ether.isSet("ldos")) {
		for (i=0;i<ether.linSize;i++) {
			aux.Ldos[i].divide(ether.iterEffective,1);
			aux.Ldos[i].print(file[4]);
		}
	}
	if (ether.isSet("chargecorrelation")) {
		vectorDivide(aux.cco,ether.iterEffective);
		for (i=0;i<4*ether.linSize;i++) {
			if (i<ether.linSize) aux.cco_aa[i%ether.linSize]=aux.cco[i];
			else if (i>=ether.linSize && i<2*ether.linSize) aux.cco_ab[i%ether.linSize]=aux.cco[i];
			else if (i>=2*ether.linSize && i<3*ether.linSize) aux.cco_ba[i%ether.linSize]=aux.cco[i];
			else aux.cco_bb[i%ether.linSize]=aux.cco[i];
		}
		vectorPrint(aux.cco_aa,"ChargeCorrelation",file[9]);
		vectorPrint(aux.cco_ab,"ChargeCorrelation",file[9]);
		vectorPrint(aux.cco_ba,"ChargeCorrelation",file[9]);
		vectorPrint(aux.cco_bb,"ChargeCorrelation",file[9]);
	}
	
	if (ether.isSet("orbitalcorrelation")) {
		vectorDivide(aux.oco,ether.iterEffective);
		for (i=0;i<4*ether.linSize;i++) {
			if (i<ether.linSize) aux.oco_aa[i%ether.linSize]=aux.oco[i];
			else if (i>=ether.linSize && i<2*ether.linSize) aux.oco_ab[i%ether.linSize]=aux.oco[i];
			else if (i>=2*ether.linSize && i<3*ether.linSize) aux.oco_ba[i%ether.linSize]=aux.oco[i];
			else aux.oco_bb[i%ether.linSize]=aux.oco[i];
		}
		vectorPrint(aux.oco_aa,"OrbitalCorrelation",file[6]);
		vectorPrint(aux.oco_ab,"OrbitalCorrelation",file[6]);
		vectorPrint(aux.oco_ba,"OrbitalCorrelation",file[6]);
		vectorPrint(aux.oco_bb,"OrbitalCorrelation",file[6]);
	}

}

template<typename MpiIoType>
void Io<MpiIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether)
{

	printSnapshot(dynVars,ether,file[2]);
}

template<typename MpiIoType>
void Io<MpiIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether,int option)
{
	string name = ether.rootname + ".last";
	std::ofstream f(name.c_str());
	printSnapshot(dynVars,ether,f);
	f.close();
}

template<typename MpiIoType>
void Io<MpiIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether,std::ofstream &f)
{
	int i,n;
	n=ether.linSize;
	
		
		
	if (ether.mpiTpemRank>0) return;			
	f<<"#Theta\n";
	for (i=0;i<n;i++)
		f<<i<<" "<<dynVars.theta[i]<<endl;

	f<<"#Phi\n";
	for (i=0;i<n;i++)
		f<<i<<" "<<dynVars.phi[i]<<endl;
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
	int j;
	for (j=0;j<ether.D;j++)
	{
		f<<"#Phonons"<<j<<"\n";
		for (i=0;i<n;i++) 
			f<<i<<" "<<dynVars.phonons[i][j]<<endl;
	}
#endif
#ifdef MODEL_KONDO_INF_TWOBANDS
	int j;
	for (j=0;j<ether.D;j++)
	{
		f<<"#Phonons"<<j<<"\n";
		for (i=0;i<n;i++) 
			f<<i<<" "<<dynVars.phonons[i][j]<<endl;
	}
#endif
#ifdef MODEL_KONDO_BCS
	unsigned int j;
	f<<"#BcsDelta\n";
	for (i=0;i<n;i++) {
		f<<i<<" "<<dynVars.bcsDelta[i]<<endl;
	}
	f<<"#BcsPhi\n";
	for (i=0;i<n;i++) {
		f<<i<<" ";
		for (j=0;j<dynVars.bcsPhi[i].size();j++) 
			f<<dynVars.bcsPhi[i][j]<<" ";
		if (ether.D==2) f<<(dynVars.bcsPhi[i][1]-dynVars.bcsPhi[i][0]);
		f<<endl;
	}
#endif				
				
}

template<typename MpiIoType>
void Io<MpiIoType>::getHostInfo()
{
				
	// Get uname info
	struct utsname *os_name;
	os_name = new utsname;
	uname(os_name);
	strcpy(hostname,os_name->sysname);
	strcat(hostname," ");
	strcat(hostname,os_name->nodename);
	strcat(hostname," ");
	strcat(hostname,os_name->release);
	strcat(hostname," ");
	strcat(hostname,os_name->version);
	strcat(hostname," ");
	strcat(hostname,os_name->machine);
	delete os_name;
}		

/*! \brief Reads the configuration file.
*/
template<typename MpiIoType>
int Io<MpiIoType>::input(char const *filename,Geometry &geometry,DynVars &dynVars,Parameters &ether,Aux &aux)
{

	int i;
	char temp2[256];
	string s,s2,s3;
	double temp=0;
	MatType tempComplex;
	int iTmp;
	vector<int> lvector;


	s=string(filename)+ ".int";
	
	if (ether.mpiRank==0) {
		if (procFile(filename,s.c_str())>0) {
			return 1;
		}
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif	
	
	
	std::ifstream fin(s.c_str());
	if (!fin || fin.bad()) {
		cerr<<"FATAL: Cannot open file: "<<filename<<endl;
		cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<endl;
		return 1;
	}
	/*! <H2>STRUCTURE OF THE INPUT FILE</H2>
	 * The configuration file is a Unix style file. 
	 * Lines beginning with # are ignored. Inline comments are not allowed.
	 * Note: Order is important here.<br>
	 * <b>Options:</b>
	
	* A single word or a comma-separated list of the following words (spaces are not allowed):
	* - \b none  You can use this tag since at least one option has to be provided. 
	*          This will not disable other options and it will simply be ignored.
	* - \b saveall  Save all classical fields after every measurement (but for for thermalization phase).
	* - \b freezephonons  In MODEL_KONDO_INF_TWOBANDS disable the evolution of phonons degrees of freedom.
	* - \b jafvector  Specify JAF as a vector of D*N doubles instead of a single number. 
	*             Allows for a spatially dependent and anisotropic JAF.
	* - \b isingspins Use ising spins instead of O(3) (heisenberg) spins.
	* - \b akw Calculate A(r,omega). With TPEM it needs postprocessing with mom2arw.pl. 
	* - \b optical Calculate sigma(omega). With TPEM it needs postprocessing. 
	* - \b havepotential Read (and use) a potential as a vector of N doubles.
	* - \b tprime Read (and use) a next-nearest neighbor hopping (only for square lattices).
	* - \b magneticfield Read (and use) a magnetic field (adds a Zeeman term to the Hamiltonian).
	* - \b groundstate Calculates the ground state and exits.
	* - \b shiftmu Shifts mu in the DOS and LDOS measurements so that the zero energy is the Fermi energy.
	* - \b cooldown Starts from a high temperature first and cools down to the value in the config. file.
	* - \b adjusttpembounds Use Lanczos to adjust TPEM bounds. (EXPERIMENTAL/UNDOCUMENTED)
	* - \b customconfig Allows for the simulation of two parallelization levels. (EXPERIMENTAL/UNDOCUMENTED)
	* - \b chargecorrelation Enables the calculation of the Quantum Charge Correlation Lattice (goes to .cco file)
	* - \b fileminimize Minimizes the number of files that are printed (caution: some information may be lost)
	* - \b spectrumbounds Allows for the specification of the spectrum bounds.
	* - \b histogrambounds Allows for the specification of histogram bounds.
	* - \b noansicolors Supresses the printing of ANSI colors to the tty.
	* - \b bcsvvector Specify bcsV as a vector of D*N doubles instead of a single number. 
	*             Allows for a spatially dependent and anisotropic BcsV.
	*/
	fin>>ether.options; //META OPTIONS none 
	/*! <b>Type of Lattice:</b> (choose a name of those between quotes)
	
	 *	- cubic types: "1d", "square", "cubic"
	 *	- all sides equal: all cubic types and "fcc" "bcc" "triangular"
	 *	- others:  "honeycomb", "prism", "rectangular"
	*/
	fin>>s; //%%INTERFACE LATTICE_TYPE
	/*! <b>Lattice Sides:</b> 
	 * A comma-separated list of numbers, one for "1d", 2 for "square",
	 * "triangular" and "honeycomb", 3 for the rest.
	 */
	fin>>s2; //%%INTERFACE LATTICE_L
	
	if (ether.isSet("verbose")) i=1;
	else i=0;
	geometry.init(s,s2,i);
	ether.linSize=geometry.volume();
	ether.D = geometry.dim();
	ether.typeoflattice = geometry.latticeName();
	
	ether.hoppings.insert(ether.hoppings.begin(),ether.D,-1);
	/*! \b Carriers:  Integer. If positive, an algorithm to adjust the chemical potential will be called to 
	 * produce the given number of electrons. If negative it is ignored. Note: (i) Even though carriers are 
	 * holes the number of electrons must be provided. (ii) In this model, 
	 * N_e + N_h = Spin_Degree_Of_Freedom*Number_Of_Bands * N. 
	 * (iii) If this option is positive then the next parameter (HAMILTONIAN_CHEMICALPOTENTIAL) 
	 * is used as the starting point to adjust the chemical potential. Ignored if TPEM is used.
	 * If non-positive it is ignored. 
	 */ 
	fin>>ether.carriers; //%%INTERFACE HAMILTONIAN_CARRIERS
	
	/*! <b>Chemical Potential: </b>
	 *	Double. Value of Chemical potential. Note: If previous parameter 
	 *	(HAMILTONIAN_NUMBER_OF_ELECTRONS > 0) then this only serves as the initial guess 
	 *	for the adjustment of the chemical potential.
	*/
	fin>>aux.varMu; //%%INTERFACE HAMILTONIAN_CHEMICALPOTENTIAL
	
	/*! \b Beta: A number that indicates the numbers to follow, followed by 
		a series of space-separated numbers, one beta per configuration
		Each is a Double. Each is the inverse of Temperature unless "temperature" option is set in 
	                 which case this must be the Temperature not its inverse.*/
	fin>>ether.numberOfBetas;
	ether.betaVector.resize(ether.numberOfBetas);
	for (size_t i = 0;i<ether.numberOfBetas;i++) {
		fin>>ether.betaVector[i]; //%%INTERFACE MONTECARLO_BETA
	
		if (ether.isSet("temperature")) {
			ether.betaVector[i] = 1.0/ether.betaVector[i];
		}
	}
	ether.beta = ether.betaVector[0];
	
	/*!  \b OutputRootname: String. Rootname for output files (See Io) */

	fin>>temp2; //%%INTERFACE OUTPUT_ROOTNAME
	ether.rootname=string(temp2); 
	
	/*! \b  MONTECARLO_THERMALIZATION: Integer. Number of Thermalization Steps */
	fin>>ether.iterTherm; //%%INTERFACE MONTECARLO_THERMALIZATION
	
	/*! \b  MONTECARLO_EFFECTIVE: Integer. Number of Measurement Steps */
	fin>>ether.iterEffective; //%%INTERFACE MONTECARLO_EFFECTIVE
	
	/*! \b  MONTECARLO_UNMEASURED: Integer. One plus the number of Steps left unmeasured during the measurement phase. */
	fin>>ether.iterUnmeasured; //%%INTERFACE MONTECARLO_UNMEASURED
	
	/*! \b MONTECARLO_WINDOW:
	 * Double. If positive or zero, this parameter determines the window used in changing the spin
	  dynamical variables (12/19/03). If negative, the value itself is ignored, and the spins are 
	  randomly changed without regard to the previous configuration. */
	fin>>ether.window; //%%INTERFACE MONTECARLO_WINDOW
	
	/*! \b MONTECARLO_FLAG: Integer. Possible values are:
	 * - \b 1 Normal Monte Carlo.
	 * - \b 0 all Monte Carlo propositions are rejected. (Useful for debugging).
	 */
	fin>>ether.mcflag; //%%INTERFACE MONTECARLO_FLAG
	
	fin>>ether.startType; //%%INTERFACE MONTECARLO_STARTTYPE
	/*! \b MONTECARLO_STARTTYPE: Integer. Starting configuration. Possible values are:
	 * - \b 0 Then theta=phi=0 for all sites. 
	 * - \b 1 Then theta and phi are chosen randomly.
	 * - \b 2 Then an antiferromagnetic configuration is used. 
	 * - \b 3 Then theta=phi=pi/2 for all sites.
	 * - \b 4 Then the program will read the DYNVARSFILENAME input parameter for the file that 
	 *        contains the initial values of the classical fields. 
	 * - \b 5 Then a CE starting config. is used (2d all sides equal lattice only).
	 * - \b 6 Like Option 4 but will also read the level at which the sav file is read.
	 */
	
	/*! \b DYNVARS_INPUT_FILE: String. The name of the input file for initial values of classical fields. 
	 * This file must be in the SAV format (see OUTPUT FILES). Provide only if STARTYPE is 4 or 6. */
	if (ether.startType==4 || ether.startType==6) fin>>ether.dynVarsInputFile; 
	
	/*! \b DYNVARS_INPUT_STARTLEVEL: Integer. The number of the set from which to read
	 * the dynvars from  DYNVARS_INPUT_FILE. Provide only if STARTYPE is 6. */
	if (ether.startType==6) {
		fin>>ether.startLevel;
	} else {
		ether.startLevel=1;
	}
	
	/*! \b HISTOGRAM_STEPS: Integer. Number of steps for the histograms. */
	fin>>ether.histSteps; //%%INTERFACE HISTOGRAM_STEPS
	
	/*! \b HISTOGRAM_BOUNDS: Two doubles. The min and max energies for the histograms.
		Provide only if option histogrambounds is set. Fine Print: 
		Applies to all histograms, in particular for  DOS and LDOS but for the optical
		conductivity  the lower bound is always zero. */
	if (ether.isSet("histogrambounds")) {
		fin>>ether.hist1;
		fin>>ether.hist2;
	} 
	
	if (ether.isSet("spectrumbounds")) {
		fin>>ether.energy1;
		fin>>ether.energy2;
	} 
	
	/*! \b HAMILTONIAN_CONCENTRATION: Integer. Number of Classical Spins. */
	fin>>ether.conc; //%%INTERFACE MONTECARLO_
	
	
	/*! \b HAMILTONIAN_BC: String. Boundary conditions. Either "periodic", "antiperiodic",
	"open" or 2D comma-separated numbers (but do not leave spaces) indicating the 
	real and imaginary part of each boundary hopping.
	*/
	fin>>s2;
	if (setBoundaryConditions(s2,ether)!=0) {
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		return 1;
	}
	
	
	/*! \b HAMILTONIAN_J: Double for  MODEL_KONDO_FINITE or two Doubles for MODEL_KONDO_DMS_MANYBANDS.
	 <i>Provide only for MODEL_KONDO_FINITE or MODEL_KONDO_DMS_MANYBANDS</i> */
	
#ifdef MODEL_KONDO_FINITE	
	fin>>temp;
	ether.J.push_back(temp);
#endif
#ifdef MODEL_KONDO_DMS_MANYBANDS
	fin>>iTmp;
	for (i=0;i<iTmp;i++) {
		fin>>temp; 
		ether.J.push_back(temp);
	}
#endif

	/*! \b HAMILTONIAN_JAF: Double or N doubles. The value of the direct exchange coupling between classical spins.
	 If option jafvector is set (see OPTIONS above) then it is a vector of Dimension*N doubles specifying 
	 Jaf[i+dir*N]. 
	 If option jafvector is not set it is a single double specifying a spatially constant value (that can be zero) for the direct
	  exchange coupling. */
	if (ether.isSet("jafvector")) {
		fin>>ether.jafFile;
		s2="#Jafvector";
		ether.JafVector.insert(ether.JafVector.begin(),ether.D*ether.linSize,0);
		loadVector(ether.JafVector,ether.jafFile,s2,1);
	} else {
		fin>>ether.numberOfJafConfigs;
		fin>>ether.jafCenter;
		fin>>ether.jafDelta;
		ether.JafVector.resize(ether.D*ether.linSize);
		//for (i=0;i<ether.D*ether.linSize;i++) ether.JafVector.push_back(temp);
	}
	/*! \b HAMILTONIAN_POTENTIAL: The name of the file containing a local potential
	 <i>but provide only if the havepotential is set (see OPTIONS above).</i>
	 */
	  i=ether.linSize;
#ifdef MODEL_KONDO_DMS_MANYBANDS
	i=2*ether.linSize; /* potentials for each band, or spin-orbit couplings first N for a, last N for b */
#endif
	ether.potential.insert(ether.potential.begin(),i,0.0);
	if (ether.isSet("havepotential")) {
		fin>>ether.potentialFile;
		s2="#Potential";
		loadVector(ether.potential,ether.potentialFile,s2,1);
	}
	
	/*! \b MAGNETIC_FIELD: 3 Doubles. A Zeeman field in Bx, By and Bz 
	 <i>but provide only if the magneticfield is set (see OPTIONS above).</i>
	 */
	 ether.magnetic.push_back(0);
	 ether.magnetic.push_back(0);
	 ether.magnetic.push_back(0);
	 
	if (ether.isSet("magneticfield")) {
		for (i=0;i<3;i++) {
			fin>>ether.magnetic[i];
		}
	}
	
	/*! \b TPRIME AND JAFPRIME: Two Doubles. A next-nearest neighbour coupling and direct exchange
	 <i>but provide only if the tprime is set (see OPTIONS above).</i>
	 */
	ether.tprime = ether.jafprime = ether.tsecond = 0.0;
	if (ether.isSet("tprime")) {
		fin>>ether.tprime; // next nearest neighbour hopping
		fin>>ether.jafprime; // next nearest neighbour direct exchange coupling
		
	}
	
	/*! \b TPEM_FLAG: Integer. Controlls the algorithm used for "diagonalization" of the one-electron sector.
	                   Possible values are:
		- 0 (diagonalization is used), 
		- 1 TPEM algorithm is used
		- 2 PEM algorithm is used
		- 3 TPEM without trace truncation is used.
	*/
	fin>>ether.tpem;
	if (ether.tpem<0 || ether.tpem>3) {
		cerr<<"Illegal option "<<ether.tpem<<" for TPEM, must be 0 (none), 1(tpem) or 2 (pem)\n";
		cerr<<"This happened here "<<__FILE__<<" "<<__LINE__<<endl;
		return 1; //ERROR CODE
	}
	
	ether.tpem_cutoff=4; // When diagonalizing this must be set also
	/*! \b TPEM_CUTOFF, TPEM_EPSPROD and TPEM_EPSTRACE: Integer, double and double respectively.
	 * The cutoff and thresholds for product and trace truncation respectively but provide these
	 * three only if TPEM_FLAG>0 */
	if (ether.tpem) {		
		fin>>ether.tpem_cutoff;
		fin>>ether.tpem_epsProd;
		fin>>ether.tpem_epsTrace;
		
	} else {
		if (ether.isSet("adjusttpembounds")) {
			cerr<<__FILE__<<": adjusttpembounds set as an option but using DIAG method\n";
			cerr<<__FILE__<<": both options are incompatible. At this point in file: "<<__LINE__<<endl;
			return 1; //ERROR CODE
		}
	}
	if (ether.carriers>0 && ether.tpem>0 && ether.mpiRank==0) {
		cerr<<"WARNING: You have selected to adjust the Chemical Potential with TPEM\n";
		cerr<<"         This is an experimental feature. This algorithm should be order N but with";
		cerr<<"         a larger pre-factor. Please report any problems.\n";
	}

	/*! \b BAND_HOPPINGS:
	 * 2x2xDIMENSION doubles. <i>Provide only for MODEL_KONDO_INF_TWOBANDS or MODEL_KONDO_DMS_MANYBANDS</i>
	 * The band hoppings in gamma, gamma', alpha, where gamma and gamma' are band indices and 
	 * alpha is x,y, etc depending on dimension. */
	 /*! \b PHONON_COUPLINGS
	  * six doubles. Provide only for MODEL_KONDO_INF_TWOBANDS The first three are the couplings 
	  * for n_i, taux_i and tauz_i and the last ones for q1_i^2, q2_i^2 and q3_i^2 as in physics reporsts
	  * 344,1. */

	
/* #ifdef MODEL_KONDO_DMS_THREEBANDS
	//cerr<<"Reading "<<9<<" band hoppings\n";
	for (i=0;i<9;i++) {
		 fin>>temp;
		 ether.bandHoppings.push_back(temp);
		// cerr<<"Read  "<<ether.bandHoppings[i]<<"\n";
	}
#endif */
	ether.numberOfOrbitals=1;	
#ifdef MODEL_KONDO_DMS_MANYBANDS
	double tmp;
	fin>>ether.numberOfOrbitals;
	iTmp=geometry.z(0);
#ifdef MODEL_KONDO_DMS_CUBIC
	iTmp=int(iTmp/2);
#endif
	cerr<<"Reading "<<ether.numberOfOrbitals<<" band hoppings\n";
	for (i=0;i<ether.numberOfOrbitals*ether.numberOfOrbitals*iTmp;i++) {
		 fin>>temp;
		 fin>>tmp;
		 tempComplex=MatType(temp,tmp);
		 ether.bandHoppings.push_back(tempComplex);
		 //cerr<<"Read  "<<ether.bandHoppings[i]<<"\n";
	}
#endif
#ifdef MODEL_KONDO_INF_TWOBANDS		
	for (i=0;i<4*ether.D;i++) {
		fin>>temp;
		tempComplex=MatType(temp,0.0);
		ether.bandHoppings.push_back(tempComplex);
	}
	for (i=0;i<3;i++) {
		fin>>temp;
		ether.phononEjt.push_back(temp);
	}
	for (i=0;i<3;i++) {
		fin>>temp;
		ether.phononEd.push_back(temp);
	}
	ether.maxPhonons=4.0/sqrt(ether.beta*maxElement(ether.phononEd));	
#endif
#ifdef MODEL_KONDO_PHONONS
	fin>>ether.phononLambda;
#endif
#ifdef MODEL_KONDO_PHONONS_EX
	fin>>ether.phononAlpha;
	fin>>ether.phononDelta;
	fin>>ether.phononGj;
#endif
	ether.bcsDelta0=0;
	ether.bcsV.insert(ether.bcsV.begin(),ether.D*ether.linSize,0);
	
#ifdef MODEL_KONDO_BCS
	fin>>ether.bcsDelta0;
	
	if (ether.isSet("bcsvvector")) {
		fin>>s3;
		s2="#BcsV";
		loadVector(ether.bcsV,s3,s2,1);
	} else {
		fin>>temp;
		for (i=0;i<ether.D*ether.linSize;i++) ether.bcsV[i]=temp;
	}	
#endif

	
	
	ether.modulus.clear();
	ether.modulus.insert(ether.modulus.begin(),ether.linSize,1);
	
	if (ether.conc>ether.linSize) {
		cerr<<"Error concentration greater than number of sites\n";
		return 1; // ERROR CODE 
	}
	
	
	/*! \b HAMILTONIAN_MODULUS: 
	 * Numbers from the set [0,N), 
	 * where N is the number of sites. These numbers are the indices of sites that have a spin.
	 * But provide only if the  
	 * following conditions apply at the same time: (i) concentration < number of sites, (ii) model is 
	 * MODEL_KONDO_FINITE and (iii) randommodulus is not set as an option. */
	if (ether.conc<ether.linSize) {
		cerr<<"conc="<<ether.conc<<" is less than linSize="<<ether.linSize<<" ==> reading modulus...\n";
		if (ether.isSet("randommodulus")) { /* It will generate a random modulus. */
			randomModulus(ether.modulus,ether.conc,ether.linSize); // from basic.h
		} else { /* It will read the file. */
			ether.modulus.assign(ether.modulus.size(),0);
			for (i=0;i<ether.conc;i++) {
				fin>>iTmp;
				if (iTmp<0 || iTmp>=ether.linSize) {
					cerr<<"Index must be >=0 and <ether.linSize but found "<<iTmp<<" instead.\n";
					cerr<<__FILE__<<" "<<__LINE__<<endl;
					return 1; //ERROR CODE
				}	
				if (!fin.good() || fin.bad()) {
					cerr<<"Premature end of input while reading spin moduli\n";
					cerr<<"This happened here "<<__FILE__<<" "<<__LINE__<<endl;
					//return 1; //ERROR CODE
				}
				ether.modulus[iTmp]=1;
			}
		}
	} // modulus reading
	
	iTmp=accumulate(ether.modulus.begin(), ether.modulus.end(),0, std::plus<int>());
	if (iTmp!=ether.conc) {
		cerr<<"Sum of spin modulus="<<iTmp<<" but concentration="<<ether.conc<<endl; 
		return 1; //ERROR CODE
	}
	if (ether.isSet("verbose")) {
		vectorPrint(ether.modulus,"modulus",cerr);
	}
	
	fin.close();

	
	return 0;
}
/*!\brief Sets the border hoppings from the information in the boundary conditions string.
	
	* \param s : can be "periodic", "antiperiodic", "open" or a comma-separated
		list of 2*dimension numbers indicating the real and imaginary
		parts of the border hoppings.
	* \param ether : the Parameters structure passed as reference. 
	*/
template<typename MpiIoType>
int Io<MpiIoType>::setBoundaryConditions(std::string const &s,Parameters &ether)
{
	int i;
	vector<double> temp;
	
	if (s=="periodic") {
		for (i=0;i<ether.D;i++) ether.hoppings[i]=MatType(-1,0);
	} else if (s=="antiperiodic") {
		for (i=0;i<ether.D;i++) ether.hoppings[i]=MatType(1,0);
	} else if (s=="open") {
		for (i=0;i<ether.D;i++) ether.hoppings[i]=MatType(0,0);
	} else {
		mysplit(s,temp,',');
		if (s.length()<2*ether.D) {
			cerr<<"Error setting boundary conditions "<<s<<endl;
			return 1;
		}
		for (i=0;i<ether.D;i++) ether.hoppings[i]=MatType(temp[2*i],temp[2*i+1]);
	} 
	return 0;
}


#endif
