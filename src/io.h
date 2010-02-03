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
#include <stdexcept>
#include "parameters.h"
#include "dynvars.h"
#include "aux.h"
#include "geometry.h"
#include "Matrix.h"
#include "SimpleReader.h"

#ifdef __LIBCATAMOUNT__ 
#include <catamount/dclock.h>
#endif

/*! \brief A class for the handling of input/output 

 * This is a singleton, that is a class that can
 * be instantiated only once. Always pass it as a reference. */
template<typename ConcurrencyIoType>
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
					if (utils::isInVector(excludedfiles,i)>=0) continue;
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
		
			template<typename SomeType>	
			void historyPrint(string const &s,const SomeType& value,int option=0)
			{
				
				SomeType result = value;
				bool debugOption = false;
				if (s=="rankGlobal=") debugOption = true;
				bool printHistory = concurrencyIo_->average(result,debugOption);
				
				if (!printHistory) return;
			
				file[0]<<s<<" "<<result;
				if (option==0) file[0]<<"\n";
			}
			
			template<typename FieldType>
			void historyPrint(string const &s,const psimag::Matrix<std::complex<FieldType> >& value,int option=0)
			{
				for (size_t p=0;p<value.n_col();p++) { // loop over pl
					file[11]<<p<<" ";
					for (size_t i=0;i<value.n_row();i++)  // loop over ks
						file[11]<<value(i,p)<<" ";
					file[11]<<"\n";
				}
			}
			
			
			void printAverages(Parameters &ether,Aux &aux);
			
			void printSnapshot(DynVars const &dynVars,Parameters const &ether);
			void printSnapshot(DynVars const &dynVars,Parameters const &ether,int option);
			void printSnapshot(DynVars const &dynVars,Parameters const &ether,std::ofstream  &f);
			void printLcd(Parameters const &ether,Aux &aux);
			void getHostInfo();
			int input(char const *filename,Geometry &geometry,DynVars &dynVars,Parameters &ether,Aux &aux);
			
			void setConcurrencyIo(ConcurrencyIoType*);
			void correctDynVarsFileName(Parameters &ether);
					
	private:
			ConcurrencyIoType* concurrencyIo_;
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


template<typename ConcurrencyIoType>
Io<ConcurrencyIoType>::Io() 
{
	
	if (isInstatiated) {
		cerr<<"FATAL: Only one Io object per program please!\n";
		cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<endl;
		exit(1);
	}
	nFiles = 12;
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
	extensions[3]=".arw_"; extensions[4]=".ldos"; extensions[5]=".opt";
	extensions[6]=".oco"; extensions[7]=".lcd"; extensions[8]=".eig";
	extensions[9]=".cco"; extensions[10]=".cor"; extensions[11]=".ncl";
			
	isInstatiated=true;
	file = new std::ofstream[nFiles];
	StartTimer();
	currentTime();
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::setConcurrencyIo(ConcurrencyIoType* ciovar)
{
	concurrencyIo_ =  ciovar;
}
		
template<typename ConcurrencyIoType>
Io<ConcurrencyIoType>::~Io()
{
	
	if (isInit) {
		writeFinalStuff();
		for (int i=0;i<nFiles;i++) {
			if (utils::isInVector(excludedfiles,i)>=0) continue;
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

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::initOutput(Parameters &ether)
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
	if (ether.mpiTpemRank==0 && concurrencyIo_->canWrite()) {			
		getHostInfo();
		getCompilerInfo();
				
		std::string filename;
		for (i=0;i<nFiles;i++) {
			//if (i>0 && rank>0) break;
			if (utils::isInVector(excludedfiles,i)>=0) continue;
			filename = std::string(ether.rootname) + extensions[i];
			file[i].open(filename.c_str());
			ether.print(file[i]);
			printMiscInfo(file[i]);
		}
		if (ether.isSet("verbose")) std::cerr<<"Io::init: success, rank="<<rank<<"\n";
		isInit=true;
	}
			
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printMiscInfo(ostream &s)
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
			
template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::writeFinalStuff()
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
		if (utils::isInVector(excludedfiles,i)>=0) continue;
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


template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::currentTime()
{
	time_t *tp = new time_t;
	time(tp);
	if (!tp) strcpy(cTime,"ERROR\n");
	else strcpy(cTime,asctime(localtime(tp)));
	int l=strlen(cTime);
	if (l>0) cTime[l-1]='\0';
	delete tp;
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printLcd(Parameters const &ether,Aux &aux)
{
	utils::vectorPrint(aux.lcd,"lcd",file[7]);
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printAverages(Parameters &ether,Aux &aux)
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
		utils::vectorDivide(aux.lcd,ether.iterEffective*ether.mpiNop2);
		utils::vectorPrint(aux.lcd,"lcd",file[7]);
		//file[0]<<"#SUMLCD="<<vectorSum(aux.lcd);
	}
	utils::vectorDivide(aux.eigM,ether.iterEffective*ether.mpiNop2);
	utils::vectorPrint(aux.eigM,"eigenvals",file[8]);
	utils::vectorDivide(aux.clasCor,ether.iterEffective*ether.mpiNop2);
	utils::vectorPrint(aux.clasCor,"clasCor",file[10]);
	if (aux.orbitalAngles.size()>0) {
		utils::vectorDivide(aux.orbitalAngles,ether.iterEffective*ether.mpiNop2);
		concurrencyIo_->vectorPrint(aux.orbitalAngles,"OrbitalAngles",file[0]);
	}
	utils::vectorDivide(aux.avMoments,ether.iterEffective*ether.mpiNop2);
	concurrencyIo_->vectorPrint(aux.avMoments,"moments",file[0]);
	
	
	for (bandIndex=0;bandIndex<ether.numberOfOrbitals;bandIndex++) {
		file[0]<<"# Band "<<bandIndex<<"\n";
		aux.Nw[bandIndex].divide(ether.iterEffective*ether.mpiNop2,1);
		aux.Nw[bandIndex].print(file[0]);
	}
	
	if (ether.isSet("optical")) {
		if (ether.tpem) {
			utils::vectorDivide(aux.opticalMoments,ether.iterEffective);
			concurrencyIo_->vectorPrint(aux.opticalMoments,"optMoments",file[5]);
		} else {
			aux.Sigma.divide(ether.iterEffective,1);
			aux.Sigma.print(file[5]);
			cerr<<"printed sigma\n";
		}
	} // optical
	
	
	if (ether.isSet("akw")) {
		if (ether.tpem) {	 
			utils::vectorDivide(aux.offdMoments,ether.iterEffective);
			utils::vectorPrint(aux.offdMoments,"moments",file[3]);
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
		utils::vectorDivide(aux.cco,ether.iterEffective);
		for (i=0;i<4*ether.linSize;i++) {
			if (i<ether.linSize) aux.cco_aa[i%ether.linSize]=aux.cco[i];
			else if (i>=ether.linSize && i<2*ether.linSize) aux.cco_ab[i%ether.linSize]=aux.cco[i];
			else if (i>=2*ether.linSize && i<3*ether.linSize) aux.cco_ba[i%ether.linSize]=aux.cco[i];
			else aux.cco_bb[i%ether.linSize]=aux.cco[i];
		}
		utils::vectorPrint(aux.cco_aa,"ChargeCorrelation",file[9]);
		utils::vectorPrint(aux.cco_ab,"ChargeCorrelation",file[9]);
		utils::vectorPrint(aux.cco_ba,"ChargeCorrelation",file[9]);
		utils::vectorPrint(aux.cco_bb,"ChargeCorrelation",file[9]);
	}
	
	if (ether.isSet("orbitalcorrelation")) {
		utils::vectorDivide(aux.oco,ether.iterEffective);
		for (i=0;i<4*ether.linSize;i++) {
			if (i<ether.linSize) aux.oco_aa[i%ether.linSize]=aux.oco[i];
			else if (i>=ether.linSize && i<2*ether.linSize) aux.oco_ab[i%ether.linSize]=aux.oco[i];
			else if (i>=2*ether.linSize && i<3*ether.linSize) aux.oco_ba[i%ether.linSize]=aux.oco[i];
			else aux.oco_bb[i%ether.linSize]=aux.oco[i];
		}
		utils::vectorPrint(aux.oco_aa,"OrbitalCorrelation",file[6]);
		utils::vectorPrint(aux.oco_ab,"OrbitalCorrelation",file[6]);
		utils::vectorPrint(aux.oco_ba,"OrbitalCorrelation",file[6]);
		utils::vectorPrint(aux.oco_bb,"OrbitalCorrelation",file[6]);
	}

}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether)
{

	printSnapshot(dynVars,ether,file[2]);
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether,int option)
{
	string name = ether.rootname + ".last";
	std::ofstream f(name.c_str());
	printSnapshot(dynVars,ether,f);
	f.close();
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::printSnapshot(DynVars const &dynVars,Parameters const &ether,std::ofstream &f)
{
	int i,n;
	n=ether.linSize;
	
		
		
	if (ether.mpiTpemRank>0) return;			
	f<<"#Theta\n";
	for (i=0;i<n;i++)
		f<<dynVars.theta[i]<<endl;

	f<<"#Phi\n";
	for (i=0;i<n;i++)
		f<<dynVars.phi[i]<<endl;
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
	int j;
	for (j=0;j<ether.D;j++)
	{
		f<<"#Phonons"<<j<<"\n";
		for (i=0;i<n;i++) 
			f<<dynVars.phonons[i][j]<<endl;
	}
#endif
#ifdef MODEL_KONDO_INF_TWOBANDS
	int j;
	for (j=0;j<ether.D;j++)
	{
		f<<"#Phonons"<<j<<"\n";
		for (i=0;i<n;i++) 
			f<<dynVars.phonons[i][j]<<endl;
	}
#endif
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::getHostInfo()
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
template<typename ConcurrencyIoType>
int Io<ConcurrencyIoType>::input(char const *filename,Geometry &geometry,DynVars &dynVars,Parameters &ether,Aux &aux)
{

	int i;
	char temp2[256];
	string s,s2,s3;
	MatType tempComplex;
	int iTmp;
	vector<int> lvector;

	typedef Dmrg::SimpleReader ReaderType;
	ConcurrencyIoType::ConcurrencyType::barrier();
	
	ReaderType reader(filename);
	
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
	* - \b nanocluster Calculates the "local structure factor" at the pre-specified momenta
	* - \b customconfig Allows for the simulation of two parallelization levels. (EXPERIMENTAL/UNDOCUMENTED)
	* - \b chargecorrelation Enables the calculation of the Quantum Charge Correlation Lattice (goes to .cco file)
	* - \b fileminimize Minimizes the number of files that are printed (caution: some information may be lost)
	* - \b spectrumbounds Allows for the specification of the spectrum bounds.
	* - \b histogrambounds Allows for the specification of histogram bounds.
	* - \b noansicolors Supresses the printing of ANSI colors to the tty.
	* - \b bcsvvector Specify bcsV as a vector of D*N doubles instead of a single number. 
	*             Allows for a spatially dependent and anisotropic BcsV.
	*/
	reader.read(ether.options); //META OPTIONS none 
	/*! <b>Type of Lattice:</b> (choose a name of those between quotes)
	
	 *	- cubic types: "1d", "square", "cubic"
	 *	- all sides equal: all cubic types and "fcc" "bcc" "triangular"
	 *	- others:  "honeycomb", "prism", "rectangular"
	*/
	reader.read(s); //%%INTERFACE LATTICE_TYPE
	/*! <b>Lattice Sides:</b> 
	 * A comma-separated list of numbers, one for "1d", 2 for "square",
	 * "triangular" and "honeycomb", 3 for the rest.
	 */
	reader.read(s2); //%%INTERFACE LATTICE_L
	
	if (ether.isSet("verbose")) i=1;
	else i=0;
	if (ether.isSet("nanocluster")) {
		reader.read(ether.plaquetteSide);
	} else {
		ether.plaquetteSide = 0;
	}
	geometry.init(s,s2,ether.plaquetteSide,i);
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
	reader.read(ether.carriers); //%%INTERFACE HAMILTONIAN_CARRIERS
	
	/*! <b>Chemical Potential: </b>
	 *	Double. Value of Chemical potential. Note: If previous parameter 
	 *	(HAMILTONIAN_NUMBER_OF_ELECTRONS > 0) then this only serves as the initial guess 
	 *	for the adjustment of the chemical potential.
	*/
	reader.read(aux.varMu); //%%INTERFACE HAMILTONIAN_CHEMICALPOTENTIAL
	
	/*! \b Beta: A number that indicates the numbers to follow, followed by 
		a series of space-separated numbers, one beta per configuration
		Each is a Double. Each is the inverse of Temperature unless "temperature" option is set in 
	                 which case this must be the Temperature not its inverse.*/
	reader.read(ether.betaVector);
	for (size_t i = 0;i<ether.betaVector.size();i++) {
		
		if (ether.isSet("temperature")) {
			ether.betaVector[i] = 1.0/ether.betaVector[i];
		}
	}
	ether.beta = ether.betaVector[0];
	
	/*!  \b OutputRootname: String. Rootname for output files (See Io) */

	reader.read(ether.rootname); 
	
	/*! \b  MONTECARLO_THERMALIZATION: Integer. Number of Thermalization Steps */
	reader.read(ether.iterTherm); //%%INTERFACE MONTECARLO_THERMALIZATION
	
	/*! \b  MONTECARLO_EFFECTIVE: Integer. Number of Measurement Steps */
	reader.read(ether.iterEffective); //%%INTERFACE MONTECARLO_EFFECTIVE
	
	/*! \b  MONTECARLO_UNMEASURED: Integer. One plus the number of Steps left unmeasured during the measurement phase. */
	reader.read(ether.iterUnmeasured); //%%INTERFACE MONTECARLO_UNMEASURED
	
	/*! \b MONTECARLO_WINDOW:
	 * Double. If positive or zero, this parameter determines the window used in changing the spin
	  dynamical variables (12/19/03). If negative, the value itself is ignored, and the spins are 
	  randomly changed without regard to the previous configuration. */
	reader.read(ether.windowVector);
	
	/*! \b MONTECARLO_FLAG: Integer. Possible values are:
	 * - \b 1 Normal Monte Carlo.
	 * - \b 0 all Monte Carlo propositions are rejected. (Useful for debugging).
	 */
	reader.read(ether.mcflag); //%%INTERFACE MONTECARLO_FLAG
	
	reader.read(ether.startType); //%%INTERFACE MONTECARLO_STARTTYPE
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
	if (ether.startType==4 || ether.startType==6) reader.read(ether.dynVarsInputFile);

	
	/*! \b DYNVARS_INPUT_STARTLEVEL: Integer. The number of the set from which to read
	 * the dynvars from  DYNVARS_INPUT_FILE. Provide only if STARTYPE is 6. */
	if (ether.startType==6) {
		reader.read(ether.startLevel);
	} else {
		ether.startLevel=1;
	}
	
	/*! \b HISTOGRAM_STEPS: Integer. Number of steps for the histograms. */
	reader.read(ether.histSteps); //%%INTERFACE HISTOGRAM_STEPS
	
	/*! \b HISTOGRAM_BOUNDS: Two doubles. The min and max energies for the histograms.
		Provide only if option histogrambounds is set. Fine Print: 
		Applies to all histograms, in particular for  DOS and LDOS but for the optical
		conductivity  the lower bound is always zero. */
	if (ether.isSet("histogrambounds")) {
		reader.read(ether.hist1);
		reader.read(ether.hist2);
	} 
	
	if (ether.isSet("spectrumbounds")) {
		reader.read(ether.energy1);
		reader.read(ether.energy2);
	} 
	
	/*! \b HAMILTONIAN_CONCENTRATION: Integer. Number of Classical Spins. */
	reader.read(ether.conc); //%%INTERFACE MONTECARLO_
	
	
	/*! \b HAMILTONIAN_BC: String. Boundary conditions. Either "periodic", "antiperiodic",
	"open" or 2D comma-separated numbers (but do not leave spaces) indicating the 
	real and imaginary part of each boundary hopping.
	*/
	reader.read(s2);
	if (setBoundaryConditions(s2,ether)!=0) {
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		return 1;
	}
	
	
	/*! \b HAMILTONIAN_J: Double for  MODEL_KONDO_FINITE or two Doubles for MODEL_KONDO_DMS_MANYBANDS.
	 <i>Provide only for MODEL_KONDO_FINITE or MODEL_KONDO_DMS_MANYBANDS</i> */
	
	reader.read(ether.J);

	/*! \b HAMILTONIAN_JAF: Double or N doubles. The value of the direct exchange coupling between classical spins.
	 If option jafvector is set (see OPTIONS above) then it is a vector of Dimension*N doubles specifying 
	 Jaf[i+dir*N]. 
	 If option jafvector is not set it is a single double specifying a spatially constant value (that can be zero) for the direct
	  exchange coupling. */
	ether.JafVector.resize(ether.linSize);
	ether.JafVector.insert(ether.JafVector.begin(),ether.D*ether.linSize,0);
	if (ether.isSet("jafdisorder")) {
		reader.read(ether.numberOfJafConfigs);
		reader.read(ether.jafCenter);
		reader.read(ether.jafDelta);
		reader.read(ether.jafSeparate);
		ether.JafVector.resize(ether.D*ether.linSize);
		//for (i=0;i<ether.D*ether.linSize;i++) ether.JafVector.push_back(temp);
	} else if (ether.isSet("jafvector")) {
		reader.read(temp2);
		ReaderType reader2(temp2);
		reader2.read(ether.JafVector); // watch of for format of external file "temp"
	} 
	
	/*! \b HAMILTONIAN_POTENTIAL: The name of the file containing a local potential
	 <i>but provide only if the havepotential is set (see OPTIONS above).</i>
	 */
	ether.potential.resize(ether.linSize);
	ether.potential.insert(ether.potential.begin(),ether.linSize,0);
	
	if (ether.isSet("potentialdisorder")) {
		reader.read(ether.numberOfMuConfigs);
		reader.read(ether.muCenter);
		reader.read(ether.muDelta);
		reader.read(ether.muSeparate);
		ether.potential.resize(ether.linSize);
	} else if (ether.isSet("havepotential")) {
		reader.read(temp2);
		ReaderType reader2(temp2);
		reader2.read(ether.potential); // watch of for format of external file "temp"
	} 
	
	/*! \b NANOCLUSTER_PARAMETERS: All integers.
	 <i>First one is the number of q's to monitor, the rest are their indices, provide only if option nanocluster option is set.</i>
	 */
	
	if (ether.isSet("nanocluster")) {
		reader.read(ether.plaquetteMeshFactor);
		reader.read(ether.q);
		ether.kmesh.init(ether.D,ether.plaquetteMeshFactor*ether.plaquetteSide);
	}
	
	/*! \b MAGNETIC_FIELD: 3 Doubles. A Zeeman field in Bx, By and Bz 
	 <i>but provide only if the magneticfield is set (see OPTIONS above).</i>
	 */
	 
	if (ether.isSet("magneticfield")) reader.read(ether.magnetic);
	
	/*! \b TPEM_FLAG: Integer. Controlls the algorithm used for "diagonalization" of the one-electron sector.
	                   Possible values are:
		- 0 (diagonalization is used), 
		- 1 TPEM algorithm is used
		- 2 PEM algorithm is used
		- 3 TPEM without trace truncation is used.
	*/
	reader.read(ether.tpem);
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
		reader.read(ether.tpem_cutoff);
		reader.read(ether.tpem_epsProd);
		reader.read(ether.tpem_epsTrace);
		
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
	
	reader.read(ether.bandHoppings); // size should be 4*ether.D	
	ether.numberOfOrbitals=1;

#ifdef MODEL_KONDO_INF_TWOBANDS
	reader.read(ether.phononEjt); // size should be 3
	reader.read(ether.phononEd); // size should be 3
	ether.maxPhonons=4.0/sqrt(ether.beta*utils::maxElement(ether.phononEd));	
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
			ether.rng.randomModulus(ether.modulus,ether.conc,ether.linSize); 
		} else { /* It will read the file. */
			reader.read(ether.modulus);
			ether.modulus.assign(ether.modulus.size(),0);
		}
	} // modulus reading
	
	iTmp=accumulate(ether.modulus.begin(), ether.modulus.end(),0, std::plus<int>());
	if (iTmp!=ether.conc) {
		cerr<<"Sum of spin modulus="<<iTmp<<" but concentration="<<ether.conc<<endl; 
		return 1; //ERROR CODE
	}
	if (ether.isSet("verbose")) {
		utils::vectorPrint(ether.modulus,"modulus",cerr);
	}

	
	return 0;
}
/*!\brief Sets the border hoppings from the information in the boundary conditions string.
	
	* \param s : can be "periodic", "antiperiodic", "open" or a comma-separated
		list of 2*dimension numbers indicating the real and imaginary
		parts of the border hoppings.
	* \param ether : the Parameters structure passed as reference. 
	*/
template<typename ConcurrencyIoType>
int Io<ConcurrencyIoType>::setBoundaryConditions(std::string const &s,Parameters &ether)
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
		utils::mysplit(s,temp,',');
		if (s.length()<size_t(2*ether.D)) {
			cerr<<"Error setting boundary conditions "<<s<<endl;
			return 1;
		}
		for (i=0;i<ether.D;i++) ether.hoppings[i]=MatType(temp[2*i],temp[2*i+1]);
	} 
	return 0;
}

template<typename ConcurrencyIoType>
void Io<ConcurrencyIoType>::correctDynVarsFileName(Parameters &ether)
{
	size_t i = ether.localRank[0];
	ether.dynVarsInputFile = ether.dynVarsInputFile+ttos(i)+".sav"; //%%INTERFACE DYNVARSINPUTROOTNAME
	ether.window = ether.windowVector[i];
}

#endif
