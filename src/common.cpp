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



#include <sys/resource.h>
#include "common.h"
#include "conductance.h"
#include "VectorGenerator.h"
#include "RandomGenerator.h"
#include "ConcurrencyParameter.h"
#include "ConcurrencyIo.h"
#if defined(CONCURRENCYMPI)
#include "ConcurrencyMpi.h"
#else
#include "ConcurrencySerial.h"
#endif
#include "Matrix.h"

template<typename ConcurrencyIoType>
bool Io<ConcurrencyIoType>::isInstatiated=false;
template<typename ConcurrencyIoType>
bool Io<ConcurrencyIoType>::isInit=false;

using namespace std;

extern void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type);
extern void createHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 MyMatrix<std::complex<double> >& matrix,Parameters const &ether,Aux &aux,int type);
extern void setupHamiltonian(MyMatrix<MatType> & matrix,Geometry const &geometry, DynVars const &dynVars, 
        Parameters const &ether,Aux &aux,int bandindex);
extern void setHilbertParams(Parameters &ether, Aux &aux, Geometry const &geometry);
extern void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry);

// template<typename ConcurrencyType>
// void setTheRankVector(Parameters& ether,std::vector<int>& v,std::vector<size_t>& w,std::vector<typename ConcurrencyType::MPIComm>& mpiCommVector)
// {
// 	
// 	size_t size0 = ether.numberOfBetas;
// 	v.resize(2);
// 	w.resize(2);
// 	w[0] = size0;
// 	w[1] = ether.mpiSize / size0;
// 	std::vector<size_t> y(2);
// 	y[0] = ether.mpiRank / w[0];
// 	y[1] = ether.mpiRank / w[1];
// 	mpiCommVector.resize(2);
// 	ConcurrencyType::MPI_Comm_split(ConcurrencyType::MPICOMMWORLD,y[0],ether.mpiRank,&(mpiCommVector[0]));
// 	ConcurrencyType::MPI_Comm_split(ConcurrencyType::MPICOMMWORLD,y[1],ether.mpiRank,&(mpiCommVector[1]));
// 	ConcurrencyType::MPI_Comm_rank(mpiCommVector[0],v[0]);
// 	ConcurrencyType::MPI_Comm_rank(mpiCommVector[1],v[1]);
// 	std::cout<<"Rank = "<<ether.mpiRank<<" v[0]="<<v[0]<<" v[1]="<<v[1]<<" size0="<<size0<<" mpiSize="<<ether.mpiSize<<"\n";
// }

/*!\brief Loads a vector from a file.
	* \param v : The vector that must have been previously allocated (output).
	* \param filename : The name of the file from where to read the vector. It should contain
		a label and then a two column index value per line with the values of
		the vector entries (input).
	* \param label : The label that indicates from where to start to read the file (input).
	* \param level : If there are many labels read the the level-th one.
	*/
int loadVector(vector<double> &v,string const &filename,string &label,int level)
{
	std::ifstream fin(filename.c_str());
	int i,n;
	string tmpLine;
	int maxline=10240;
	char *buffer;
	
	buffer  = new char[maxline];
		
	if (!fin || fin.bad()) {
		cerr<<"Cannot open file "<<filename<<endl;
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1); // throw an exception instead, caller must catch it
	}
	
	//cerr<<"load vector: I just successfully opened: "<<filename<<endl;
	for (i=0;i<level;i++) { // .sav will have 2 data sets, read the level-th data set
		while(!fin.eof()) {
			fin.getline(buffer,maxline);
			tmpLine=string(buffer);
			if (fin.fail()) {
				cerr<<"loadVector: Premature end of file="<<filename<<" no label="<<label<<" found.\n";
				cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
				exit(1);
			}
			//cerr<<"Load vector, before chopping "<<tmpLine<<"*"<<endl;
			utils::mychop(tmpLine);
			//cerr<<"Load vector, after chopping "<<tmpLine<<"*"<<endl;
			if (tmpLine==label) {
				//cerr<<"I read "<<tmpLine<<" which IS equal to "<<label<<endl;
				break;
			}
			//cerr<<"I read "<<tmpLine<<" which is not equal to "<<label<<endl;
		}
	}
	
	if (!(tmpLine==label)) {
		cerr<<"loadVector: Problem trying to read file "<<filename<<endl;
		cerr<<"Label="<<label<<" not found\n";
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	n = v.size();
	vector<double> c;
	
	for (i=0;i<n;i++) {
		fin.getline(buffer,maxline);
		tmpLine=string(buffer);
		utils::mysplit(tmpLine,c,' ');
		if (i!=c[0]) {
			cerr<<"loadVector: Mismatch while reading label="<<label<<" in file "<<filename;
			cerr<<" "<<i<<" not equal "<<c[0]<<endl;
			cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
			exit(1);// throw an exception instead, caller must catch it
		} 
		v[i]=c[1];
	}
	fin.close();
	delete [] buffer;
	return 0;
}


void tmpValues(double &a,double &b,double &mu,double &beta,int option)
{
	static double a1,b1,mu1,beta1;
	if (option==0) { // set static members
		a1 = a;
		b1 = b;
		mu1 = mu;
		beta1 = beta;
	} else { //retrieve static members
		a = a1;
		b = b1;
		mu = mu1;
		beta = beta1;
	}
}
	

// should go elsewhere:
template<typename T>
T volumeOf(std::vector<T>& v)
{
	if (v.size()==0) return 0;
	T tmp = v[0];
	for (size_t i=1;i<v.size();i++) tmp *= v[i];
	return tmp;
}

void setTheSizeVector(Parameters& ether,std::vector<size_t>& w)
{
	
	w.resize(3);
	w[0] = ether.betaVector.size();
	w[1] = ether.numberOfJafConfigs;
	w[2] = ether.numberOfMuConfigs;
	
	if (size_t(ether.mpiSize)!=volumeOf(w)) throw std::runtime_error("Number of processors mismatch\n");
}

template<typename ConcurrencyType>
void setTheRankVector(Parameters& ether,std::vector<int>& ranks,const std::vector<size_t>& sizes,std::vector<typename ConcurrencyType::MPIComm>& mpiCommVector)
{
	size_t n = sizes.size();
	ranks.resize(n);
	std::vector<size_t> keysForSplit(n);
	mpiCommVector.resize(n);
	keysForSplit[0] = ether.mpiSize;
	
	std::vector<size_t> volume(n+1);
	volume[n] = 1;
	for (int i=n - 1;i>=0;i--) 
		volume[i] = volume[i+1]*sizes[i];
		 
	for (size_t i=0;i<n;i++) { 
		if (i<n-1) keysForSplit[i] = ether.mpiRank % volume[i+1] + size_t(ether.mpiRank/volume[i])*sizes[i+1];
		else keysForSplit[i] =  ether.mpiRank % volume[i+1] + size_t(ether.mpiRank / volume[i]);
		ConcurrencyType::MPI_Comm_split(ConcurrencyType::MPICOMMWORLD,keysForSplit[i],ether.mpiRank,&(mpiCommVector[i]));
		ConcurrencyType::MPI_Comm_rank(mpiCommVector[i],ranks[i]);
		std::cout<<"\t ranks["<<i<<"]="<<ranks[i]<<" size["<<i<<"]="<<sizes[i]<<"Rank = "<<ether.mpiRank<<" mpiSize="<<ether.mpiSize<<"\n";;
	}
	
}

template<typename ConcurrencyType>
void registerHook(Parameters& ether,ConcurrencyIo<ConcurrencyType>** ciovar)
{
	std::vector<typename ConcurrencyType::MPIComm> mpiCommVector;
	setTheSizeVector(ether,ether.localSize);
	setTheRankVector<ConcurrencyType>(ether,ether.localRank,ether.localSize,mpiCommVector);
	size_t counter = 0;
	
	//example of non-random
	VectorGenerator<double> betaGenerator(ether.betaVector,ether.localRank[counter]);
	//VectorGenerator<string> dynVarsFileGenerator(ether.dynVarsInputFileVector,ether.localRank[counter]);
	typedef ConcurrencyParameter<double,VectorGenerator<double>,Parameters,ConcurrencyType> ConcurrencyParameterBeta;
	//typedef ConcurrencyParameter<string,VectorGenerator<string>,Parameters,ConcurrencyType> ConcurrencyParameterDynvarsFile;
	ConcurrencyParameterBeta beta(ether.beta,ether,betaGenerator,ConcurrencyParameterBeta::SEPARATE,
				      ether.localRank[counter],mpiCommVector[counter],ether.localSize[counter]); // beta 4 10 20 30 40
	//ConcurrencyParameterDynvarsFile dynVarsInputFile(ether.dynVarsInputFile,ether,dynVarsFileGenerator,ConcurrencyParameterDynvarsFile::GATHER,
	//						 ether.localRank[counter],mpiCommVector[counter],ether.localSize[counter]);
	counter++;
	
	// example of random
	typedef RandomGenerator<double> RandomGeneratorType;
	RandomGeneratorType rng("bimodal");
	
	typedef ConcurrencyParameter<std::vector<double>,RandomGeneratorType,Parameters,ConcurrencyType> ConcurrencyParameterRandType;
	rng.setBimodal(ether.jafCenter,ether.jafDelta); 
	rng.seed(ether.localRank[counter]);
	ConcurrencyParameterRandType jafvector(ether.JafVector,ether,rng,ether.jafSeparate,
					  ether.localRank[counter],mpiCommVector[counter],ether.localSize[counter]); // uses the sequence
	counter++;
	
	// example of random
	rng.setBimodal(ether.muCenter,ether.muDelta); 
	rng.seed(ether.localRank[counter]);
	ConcurrencyParameterRandType muvector(ether.potential,ether,rng,ether.muSeparate,
					  ether.localRank[counter],mpiCommVector[counter],ether.localSize[counter]); //uses the sequence

	
	// allocate
	ciovar[0] = new ConcurrencyIo<ConcurrencyType>(beta,jafvector,muvector);
}

int spf_entry(int argc,char *argv[],int mpiRank=0, int mpiSize=1)
{
	Parameters ether;
	char sfile[256];
	TpemOptions tpemOptions;
	
	tpemOptions.mpi_nop2 =1;
	
	rlimit myrlim;
        myrlim.rlim_cur=myrlim.rlim_max=0;
        setrlimit(RLIMIT_CORE, &myrlim);

        tpem_init(argc,argv,tpemOptions,sfile);
	ether.mpiRank =mpiRank;
	ether.mpiNop1 = tpemOptions.mpi_nop1;
	ether.mpiNop2 = tpemOptions.mpi_nop2;
	ether.mpiTpemRank = tpemOptions.tpem_rank;
	ether.mpiSize = mpiSize;
	
	//cerr<<"main: global rank="<<tpemOptions.rank<<" local rank="<<tpemOptions.tpem_rank;
	//cerr<<" nop1="<<ether.mpiNop1<<" "<<ether.mpiNop2<<endl;
	

	if (argc<2) {
		cerr<<"Not enough arguments. USAGE IS\n";
		cerr<<argv[0]<<" filename\n";
		return 1;
	}
	
	int iter,iter2;
	Geometry geometry;
	DynVars dynVars;
	Aux aux;
#if defined(CONCURRENCYMPI)
	typedef ConcurrencyIo<ConcurrencyMpi> ConcurrencyIoType;
#else
	typedef ConcurrencyIo<ConcurrencySerial> ConcurrencyIoType;
#endif
	
	Io<ConcurrencyIoType> io;
	
	srand(time(0));

	
	if (io.input(sfile,geometry,dynVars,ether,aux)) {
		if (ether.mpiRank==0) cerr<<"A problem ocurred while trying to load file "<<sfile<<endl;
		return 1;
	}
	
	ConcurrencyIoType** ConcurrencyIo = new ConcurrencyIoType*[1];
	registerHook(ether,ConcurrencyIo);
	io.setConcurrencyIo(ConcurrencyIo[0]);

	// has to happend after registerHook because we use ether.localRank that is setup by registerHook
	// and            before setupVariables that uses it.
	io.correctDynVarsFileName(ether);
	
	// enable custom config
	if (ether.mpiNop2>1 && (ether.isSet("optical") || ether.isSet("akw"))) {
		if (ether.mpiRank==0) {
			cerr<<"main: optical and akw is set and conflicts with grid computing\n";
			cerr<<"main: tip: Either remove the options optical and/or akw or disable grid computing ";
			cerr<<" by not providing a second argument in the command line\n";
		}
		tpem_finalize();
		return 1;
	}
	if (ether.mpiNop2>1) {
		double shift=0.0;
		if (argc>=5) shift=atof(argv[4]);
		customConfig(ether,aux,tpemOptions,atoi(argv[3]),shift);
	}
	
	setupVariables(geometry,dynVars,ether,aux);
	tpemOptionsFill(tpemOptions,ether);
	
	io.initOutput(ether);

	setupHamiltonian(aux.matrix,geometry,dynVars,ether,aux,0);
	
	if (ether.tpem>0 && ether.isSet("sineupdate")) {
		cerr<<"You have selected the TPEM and sineupdate.\n";
		cerr<<"This is not currently supported.\n";
		cerr<<"At this point: "<<__FILE__<<" : "<<__LINE__<<endl;
		return 1;
	} 
	  if (ether.isSet("groundstate")) {
                calcGroundState(geometry,dynVars,ether,aux);
                io.printSnapshot(dynVars,ether);
                return 0;
        }
	if (ether.mpiRank==0) {	
		ether.welcome();
		io.fileOut()<<"\nTHERMALIZATION PHASE ----------\n";
		io.fileOut()<<"PROGRESS:            "<<endl;
		cout.flush();
	}

	if (ether.isSet("cooldown")) {
		double betaMax = ether.beta;
		int therm = ether.iterTherm;

		for (iter=0;iter<therm;iter++) {
			ether.beta = (iter+1.0)*betaMax/therm;  
			doMonteCarlo(geometry,dynVars,ether,aux,tpemOptions);
		}
	}
	
	  
	for (iter=0;iter<ether.iterTherm;iter++) {
		utils::printProgress(iter,ether.iterTherm,10,'*',ether.mpiRank);
		doMonteCarlo(geometry,dynVars,ether,aux,tpemOptions);
	}
	if (ether.mpiRank==0) {
		io.fileOut()<<"\nMEASUREMENT PHASE ----------\n";
		io.fileOut()<<"PROGRESS:            "<<endl;
		cout.flush();
	}

	for (iter=0;iter<ether.iterEffective;iter++) {
		utils::printProgress(iter,ether.iterEffective,10,'*',ether.mpiRank);

		
		for (iter2=0;iter2<ether.iterUnmeasured;iter2++) {
			doMonteCarlo(geometry,dynVars,ether,aux,tpemOptions);
		}
		doMeasurements(iter,dynVars,geometry,io,ether,aux,tpemOptions);
	}
	
	io.printAverages(ether,aux);
	io.printSnapshot(dynVars,ether);
	io.fileOut()<<"\n\n";
	io.fileError()<<"Done!\a\n";	
	tpem_finalize();
	

	return 0;
}

void tpem_matrix_print (tpem_sparse *t)
{
        unsigned int i, k;

        for (i = 0; i < t->rank; i++) {
                for (k = t->rowptr[i]; k < t->rowptr[i + 1]; k++) {
                        cerr<<"("<<real(t->values[k])<<","<<imag(t->values[k])<<") ";
                }
		cerr<<endl;
        }
}

void tpemOptionsFill(TpemOptions &tpemOptions,Parameters const &ether)
{
	
	tpemOptions.cutoff = ether.tpem_cutoff;
	if (ether.tpem==1) {
		tpemOptions.tpemType = "tpem";
	} else if (ether.tpem==2) {
		tpemOptions.tpemType ="pem";
	}
	tpemOptions.coeffs="fast";
#ifndef NO_GSL		
	tpemOptions.coeffs = "gsl";
#endif
	if (ether.carriers>0) tpemOptions.coeffs="fast";
	tpemOptions.epsTrace = ether.tpem_epsTrace;
	tpemOptions.epsProd = ether.tpem_epsProd;
}



#ifdef MODEL_KONDO_PHONONS
double electronPhononTerm(int p,Geometry const &geometry, DynVars const &dynVars,Parameters const &ether)
{
	double sum=0;
	int alpha,j;
	for (alpha=0;alpha<geometry.z(p);alpha++) {
		if (alpha>=geometry.dim()) break;
		j = geometry.neighbor(p,2*alpha+1);
		sum += dynVars.phonons[j][alpha];
		sum -= dynVars.phonons[p][alpha];
	}
	return sum;
}

#endif

void customConfigAlter(int x, int y,Parameters &ether,Aux &aux)
{
	ether.beta = 1.0/ether.beta;
	ether.beta = ether.beta * (1+y);
	ether.beta = 1.0/ether.beta;
	if (ether.isSet("ccdisorder")) {
                // use x for disorder
                ether.jafFile = ether.jafFile + utils::ttos(y) + ".txt";
                string s2="#Jafvector";
                loadVector(ether.JafVector,ether.jafFile,s2,1);
                ether.options = ether.options + ",jafvector";
        }
	
}




void customConfig(Parameters &ether,Aux &aux,TpemOptions const &tpemOptions,
int nop3,double shift)
{

        if (!ether.isSet("customconfig")) {
                if (ether.mpiRank==0) {
                        cerr<<"customConfig: I was expecting customconfig option to be set.\n";
                        cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
                }
                exit(1);
        }
        // change an input property
        int color = int(tpemOptions.rank/tpemOptions.mpi_nop1);
        int x,y;

        if (tpemOptions.mpi_nop2< nop3 || tpemOptions.mpi_nop2 % nop3!=0) {
                if (ether.mpiRank==0) {
                        cerr<<"n_params does not divide tpemOptions.mpi_nop2\n";
                        cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
                }
                exit(1);
        }

 
        y=int(color/nop3);
        x=color % nop3;
	
	customConfigAlter(x,y,ether,aux);
	
        // change rootname
        std::stringstream ss;
        std::string str;
        ss << color;
        ss >> str;
        ether.rootname = ether.rootname + str;

}


double Zeeman(DynVars const &dynVars,Geometry const &geometry,Parameters const &ether)
{
	vector<double> s(3,0);
	int i;
	double ret=0;
	
	s[0]=s[1]=s[2]=0;
	for (i=0;i<geometry.volume();i++) {
		s[0]+= ether.modulus[i]*sin(dynVars.theta[i])*cos(dynVars.phi[i]);
		s[1]+= ether.modulus[i]*sin(dynVars.theta[i])*sin(dynVars.phi[i]);
		s[2]+= ether.modulus[i]*cos(dynVars.theta[i]);
	}
	for (i=0;i<3;i++) ret += s[i] * ether.magnetic[i];
	
	return ret;
}

	
double dSDirect(DynVars const &dynVars, DynVars const &dynVars2, int i, Geometry const &geometry, Parameters const &ether)
{
	int j,k;
	double dS=0;
	double tmp,jaf;
	
	if (!geometry.isCubicType()) return 0;
		
	for (k = 0; k<geometry.z(i); k++){
		j=geometry.neighbor(i,k);
		tmp = (sin(dynVars2.theta[i])*cos(dynVars2.phi[i])-sin(dynVars.theta[i])*cos(dynVars.phi[i]))*sin(dynVars.theta[j])*cos(dynVars.phi[j]) 
		    + (sin(dynVars2.theta[i])*sin(dynVars2.phi[i])-sin(dynVars.theta[i])*sin(dynVars.phi[i]))*sin(dynVars.theta[j])*sin(dynVars.phi[j]) 
		    + (cos(dynVars2.theta[i])-cos(dynVars.theta[i]))*cos(dynVars.theta[j]);
		if (k==0 || k%2==0) jaf=ether.JafVector[i+k*ether.linSize/2];
		else jaf=ether.JafVector[j+(k-1)*ether.linSize/2];
		
		dS += jaf*tmp;
	}
	if (ether.isSet("magneticfield")) tmp = Zeeman(dynVars2,geometry,ether)-Zeeman(dynVars,geometry,ether);
	else tmp =0;
	
	dS += tmp;
	

	return dS;
}




double dPhonons(DynVars const &dynVars,DynVars const &dynVars2,int i,Geometry const &geometry, Parameters const &ether,const Phonons<Parameters,Geometry>& phonons)
{
	double dE=0.0;
	
#ifdef MODEL_KONDO_INF_TWOBANDS
	for (size_t alpha=0;alpha<dynVars.phonons[i].size();alpha++) {
		double tmp=0;
		tmp += square(phonons.calcPhonon(i,dynVars2,alpha));
		tmp -= square(phonons.calcPhonon(i,dynVars,alpha));
		for (int k=0;k<geometry.z(i);k++) {
			int j = geometry.neighbor(i,k);
			tmp += square(phonons.calcPhonon(j,dynVars2,alpha));
			tmp -= square(phonons.calcPhonon(j,dynVars,alpha));
		}
		dE += ether.phononEd[alpha]*tmp;
	}
#else	
	for (size_t j=0;j<dynVars.phonons[i].size();j++)
		dE += square(dynVars2.phonons[i][j])-square(dynVars.phonons[i][j]);
#endif
	return dE;
}	

double directExchange2(DynVars const &dynVars, Geometry const &geometry, Parameters const &ether)
{
	int i,j,k;
	double dS=0;
	double tmp;
	int n = geometry.volume();
	double t1,t2,p1,p2;
	
	for (i=0;i<n;i++) {
	for (k = 0; k<geometry.z(i,2); k++){
		j=geometry.neighbor(i,k,2); /**next nearest neighbor */
		t1=dynVars.theta[i];
		t2=dynVars.theta[j];
		p1=dynVars.phi[i];
		p2=dynVars.phi[j];
		tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
		//dS += ether.jafprime*tmp; // FIXME
	}}
	
	return dS*0.5;
}

double calcSuperExchange(DynVars const &dynVars,Geometry const &geometry, Parameters const &ether)
{
	int i,j,k;
	double dS=0;
	double jaf,tmp;
	int n = geometry.volume();
	double t1,t2,p1,p2;
	
	if (!geometry.isCubicType()) return 0;
	
	for (i=0;i<n;i++) {
	for (k = 0; k<geometry.z(i); k++){
		j=geometry.neighbor(i,k);
		t1=dynVars.theta[i];
		t2=dynVars.theta[j];
		p1=dynVars.phi[i];
		p2=dynVars.phi[j];
		if (k==0 || k%2==0) jaf=ether.JafVector[i+k*ether.linSize/2];
		else jaf=ether.JafVector[j+(k-1)*ether.linSize/2];
		tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
		dS += jaf*tmp;
	}}
	
	return dS*0.5;
}

double calcPhononEnergy(DynVars const &dynVars, Geometry const &geometry, const Parameters& ether,Phonons<Parameters,Geometry> const &phonons)
  // Written by IS Oct-19-04.
  // Calculates pure phononic energy of the system (last term in the Hamiltonian).
{
  size_t n = geometry.volume();
  double dE=0.0;
	
  for (size_t i=0;i<n;i++) 
    {

#ifdef MODEL_KONDO_INF_TWOBANDS
      for (size_t alpha=0;alpha<dynVars.phonons[i].size();alpha++) 
	{
	  double tmp = phonons.calcPhonon(i,dynVars,alpha);
	  // cerr << "site "<< i << "phonon "<< alpha << ": " <<tmp << endl;
	  dE += ether.phononEd[alpha]*square(tmp);
	}
#else	
      for (size_t j=0;j<dynVars.phonons[i].size();j++)
	dE += square(dynVars.phonons[i][j]);
#endif
    }
  return dE; //no 0.5 factor here!
}

#ifdef MODEL_KONDO_PHONONS

double calcPhononV(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether)
{
	int i;
	double sum=0;
	for (i=0;i<ether.linSize;i++) {
		sum+= electronPhononTerm(i,geometry,dynVars,ether);	
	}
	return sum;
}

double calcPhononV2(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether)
{
	int i;
	double sum=0;
	for (i=0;i<ether.linSize;i++) {
		sum += square(electronPhononTerm(i,geometry,dynVars,ether));	
	}
	return sum;
}

#endif	

void kTpemMoments(vector<double> const &moments,Aux &aux, Parameters const &ether)
{
	
	for (unsigned int m=0;m<moments.size();m++) {
		aux.avMoments[m] += moments[m];
	}
}	



void crsToFull (MyMatrix<MatType> &a,tpem_sparse const *matrix,Parameters const &ether,Aux &aux)
{
	unsigned int i,n,k;
	
	n=a.getRank();
	
	for (i = 0; i < n ; i++) for (k=0;k<n;k++) a.set(i,k,0.0);
		
	for (i = 0; i < n; i++){
		for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++){			
			a.set(matrix->colind[k],i,matrix->values[k]);
		}
	}
	//if (!a.isHermitian()) {
	//	std::cerr<<a<<endl;
	//	exit(1);
	//}
	
}

double calculate_dS (vector<double> &moment,tpem_sparse *matrix1, tpem_sparse *matrix, 
	vector<unsigned int> const &support,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff = ether.tpem_cutoff;
        static vector<double>  coeffs(cutoff),moment0(cutoff),moment1(cutoff);
        double dS;
        static int firstcall=1;
        double a=aux.varTpem_a,b=aux.varTpem_b,mu=aux.varMu,beta=ether.beta;
	double e1=aux.varTpem_b-aux.varTpem_a, e2=aux.varTpem_a+aux.varTpem_b;
	double e11=aux.btmp-aux.atmp, e21=aux.atmp+aux.btmp;
	
	if (e1 > e11) e1 = e11;
	if (e2 < e21) e2 = e21;
	
	
        if (firstcall || ether.carriers>0) {
		firstcall=0;
                tmpValues(a,b,mu,beta,0);
		tpem_calculate_coeffs (coeffs,effective_action,tpemOptions);
        } 
       	if (ether.isSet("adjusttpembounds")) {
		a=0.5*(e2-e1);
		b=0.5*(e2+e1);
		tpem_calculate_coeffs (coeffs,effective_action,tpemOptions);
	}
	
		
	tmpValues(a,b,mu,beta,0);
        
	tpem_sparse *mod_matrix1, *mod_matrix;
	
	if (ether.isSet("adjusttpembounds")) {
		// full copy/allocation
		mod_matrix1 = new_tpem_sparse(ether.hilbertSize,matrix1->rowptr[ether.hilbertSize]);
		tpem_sparse_copy (matrix1, mod_matrix1);
		mod_matrix = new_tpem_sparse(ether.hilbertSize,matrix->rowptr[ether.hilbertSize]);
		tpem_sparse_copy (matrix, mod_matrix);
	
		tpem_sparse_scale(mod_matrix1,a/aux.atmp,(b-aux.btmp)*aux.atmp/(a*a));
		tpem_sparse_scale(mod_matrix,a/aux.varTpem_a,(b-aux.varTpem_b)*aux.varTpem_a/(a*a));
		
	} else {		
		// just pointers
		mod_matrix1 = matrix1;
		mod_matrix = matrix;
	}
	
	tpem_calculate_moment_diff (mod_matrix1, mod_matrix,moment,support, tpemOptions);
	
        dS = -tpem_expansion (coeffs, moment);
		 	
	if (ether.isSet("adjusttpembounds")) {
		tpem_sparse_free(mod_matrix);
		tpem_sparse_free(mod_matrix1);
	}
	
	
	return dS;
}




/***************** DIAGONALIZATION WRAPPERS ************************/
void diag(vector<double> &eig,Geometry const &geometry,DynVars const &dynVars,
	Parameters const &ether,Aux &aux,char jobz)
{
	if (eig.size() != size_t(ether.hilbertSize)) {
		cerr<<"Internal Error at diag\n";
		exit(1);
	}
	if (jobz=='v') jobz='V';
	
	setupHamiltonian(aux.matrix,geometry,dynVars,ether,aux,0);
	if (ether.isSet("matrixprint")) {
		cerr<<aux.matrix<<endl;
		//matrixPrint(matrix);
	}

	if (ether.isSet("blockmatrix")) {
		// diagonalizes the matrix assuming it is block diagonal
		MyMatrix<MatType> lowerMatrix,upperMatrix;
		size_t blocksize=aux.matrix.getRank()/2;
		vector<double> eig1(blocksize),eig2(blocksize);
		
		
		lowerMatrix.init(blocksize,0);
		upperMatrix.init(blocksize,0);
		for (size_t i=0;i<blocksize;i++) {
			for (size_t j=0;j<blocksize;j++)  {
				upperMatrix(i,j)=aux.matrix(i,j);
				lowerMatrix(i,j)=aux.matrix(i+blocksize,j+blocksize);
			}
		}
		diag(upperMatrix,eig1,jobz);
		diag(lowerMatrix,eig2,jobz);
		for (size_t i=0;i<blocksize;i++) {
			eig[i]=eig1[i];
			eig[i+blocksize]=eig2[i];
		}
		if (jobz=='V') { // eigenvectors
			for (size_t i=0;i<blocksize;i++) {
				for (size_t j=0;j<blocksize;j++) {
					aux.matrix(i,j)=upperMatrix(i,j);
					aux.matrix(i,j+blocksize)=0;
					aux.matrix(i+blocksize,j)=0;
					aux.matrix(i+blocksize,j+blocksize)=lowerMatrix(i,j);
				
				}
			}
		}
		
		
	} else {
		diag(aux.matrix,eig,jobz);
	}
	
	if (jobz!='V') sort(eig.begin(), eig.end(), less<double>());
	if (ether.isSet("eigprint")) {
		for (size_t i=0;i<eig.size();i++) {
			cerr<<"eig["<<i<<"]="<<eig[i]<<endl;
		}
		exit(1);
	}
}


/***************** TPEM SPECIFIC CODE********************************/
bool kTpemAllBands(int i,Geometry const &geometry,DynVars const &dynVars,
	Parameters const &ether,Aux &aux,double dsDirect,TpemOptions const &tpemOptions)
{
	
	
	double dS=0;
	bool accept;
	vector<unsigned int> support;
	bool glauber=true;
	vector<double> moment(ether.tpem_cutoff);

	dS =0;
	setSupport(support,i,geometry);
	
	kTpemHamiltonian (geometry, dynVars,aux.sparseTmp[0],ether,aux,0);
	if (ether.isSet("adjusttpembounds")) {
		if (tpemAdjustBounds(aux.sparseTmp[0],ether,aux)!=0) {
			cerr<<"Cannot adjust bounds for tpem spectrum\n";
			exit(1);
		}
	}
	dS += calculate_dS (moment,aux.sparseMatrix[0], aux.sparseTmp[0],
	support,ether,aux,tpemOptions);

	dS -= ether.beta*dsDirect;
	tpem_bcast(&dS,tpemOptions);
	
	if (glauber) {
		if (dS<0) {
			dS=exp(dS)/(1.0+exp(dS));
		} else {
			dS=1.0/(1.0+exp(-dS));
		}
		if (dS>ether.rng.myRandom()) accept=true;
		else accept=false;
	}	else {
		accept = (dS > 0.0 || ether.rng.myRandom () < exp (dS));
	}
	if (accept) {
		// update current moments
		for (int j=0;j<ether.tpem_cutoff;j++) aux.curMoments[j] -= moment[j];		
	}
	return accept;
}

void kTpemSetAllBands(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,Aux &aux)
{
	
	
	
	kTpemHamiltonian (geometry, dynVars,aux.sparseMatrix[0],ether,aux,0);
	if (ether.isSet("adjusttpembounds")) {
		if (tpemAdjustBounds(aux.sparseMatrix[0],ether,aux)!=0) {
			cerr<<"Cannot adjust bounds for tpem spectrum\n";
			exit(1);
		}
	}
}


// TPEM_FIXME: vector allocations are for moments, matrix allocation 
//                                                 should be for sparse matrix	
void setupVariables(Geometry const &geometry,DynVars &dynVars,Parameters &ether,Aux &aux)
{
	int i,d,n,matsize,j;
	n=geometry.volume();
	d=geometry.dim();
	vector<double> vtmp;
	
	if (ether.isSet("randomize")) ether.rng.myRandomSeed(time(0));
	
	// NOTE: THIS FUNCTION SHOULD NOT BE REENTRANT--> FIXME
		
	setHilbertParams(ether,aux,geometry);
	
	aux.nac.insert(aux.nac.begin(),ether.classFieldList.size(),0);
	
	dynVars.theta.insert(dynVars.theta.begin(),n,0.0);
	dynVars.phi.insert(dynVars.phi.begin(),n,0.0);
	
	dynVars.Tx.insert(dynVars.Tx.begin(),n,0.0);
	dynVars.Ty.insert(dynVars.Ty.begin(),n,0.0);
	dynVars.Tz.insert(dynVars.Tz.begin(),n,0.0);
	
	for (i=0;i<d;i++) vtmp.push_back(0.0);
	
	
	
	for (i=0;i<n;i++) {
		dynVars.phonons.push_back(vtmp);
		
	}
	for (i=0;i<n;i++) {
		switch (ether.startType) {
			case 1:
				dynVars.theta[i]=M_PI*ether.rng.myRandom();
				dynVars.phi[i]=2*M_PI*ether.rng.myRandom();
				
				if (!ether.isSet("freezephonons")) {
					if (ether.isSet("verbose") && ether.mpiRank==0) cerr<<"Randomizing phonons\n";
					for (j=0;j<d;j++) dynVars.phonons[i][j]=ether.maxPhonons*(ether.rng.myRandom()-0.5);
				}								
				break;
			case 2:
				dynVars.theta[i]=(utils::parity(i,d,geometry.length())==1) ? 0 : M_PI;
				break;
			case 3:
				dynVars.theta[i]=M_PI*0.5;
				dynVars.phi[i]=M_PI*0.5;
				break;
		}
		if (ether.isSet("isingspins")) {
			if (dynVars.theta[i]>M_PI/2) dynVars.theta[i]=M_PI;
			else dynVars.theta[i]=0;
			dynVars.phi[i]=0;
		}
	}
	
	aux.eigM.insert(aux.eigM.begin(),ether.hilbertSize,0.0);
	aux.lcd.insert(aux.lcd.begin(),n,0.0);
	aux.clasCor.insert(aux.clasCor.begin(),n,0.0);

	aux.cco_aa.insert(aux.cco_aa.begin(),n,0.0);
	aux.cco_ab.insert(aux.cco_ab.begin(),n,0.0);
	aux.cco_ba.insert(aux.cco_ba.begin(),n,0.0);
	aux.cco_bb.insert(aux.cco_bb.begin(),n,0.0);
	aux.cco.insert(aux.cco.begin(),n*4,0.0);
	aux.oco_aa.insert(aux.oco_aa.begin(),n,0.0);
	aux.oco_ab.insert(aux.oco_ab.begin(),n,0.0);
	aux.oco_ba.insert(aux.oco_ba.begin(),n,0.0);
	aux.oco_bb.insert(aux.oco_bb.begin(),n,0.0);
	aux.oco.insert(aux.oco.begin(),n*4,0.0);
	
	aux.bcsCorxx.insert(aux.bcsCorxx.begin(),n,0.0);
	aux.bcsCorxy.insert(aux.bcsCorxy.begin(),n,0.0);
	if (ether.isSet("orbitalangles")) aux.orbitalAngles.insert(aux.orbitalAngles.begin(),n,0.0);
		
	// Allocate the matrix
	matsize=ether.hilbertSize;
   	if (ether.tpem) {
		aux.sparseMatrix = new tpem_sparse*[1];
		aux.sparseTmp = new tpem_sparse*[1];
		aux.sparseMatrix[0] = new_tpem_sparse (matsize,ether.nonzero);
		aux.sparseTmp[0] = new_tpem_sparse (matsize, ether.nonzero);
		 if (ether.isSet("conductance")) aux.matrix.init(matsize,0.0);
	} else {
		aux.matrix.init(matsize,0.0);
	}
	if (!ether.isSet("histogrambounds")) {
		ether.hist1=ether.energy1;
		ether.hist2=ether.energy2;
	}
	
	// FIXME: When bands are shifted remember to change normalization
	
	aux.eigOneBand.insert(aux.eigOneBand.begin(),matsize,0.0);
	aux.eigAllBands.insert(aux.eigAllBands.begin(),matsize,0.0);
	
	for (i=0;i<ether.numberOfOrbitals;i++) aux.Nw[i].init(ether.hist1,ether.hist2,ether.histSteps);
	
	
	if (ether.isSet("akw")) {
		j=n;
		if (ether.typeofmodel=="MODEL_KONDO_DMS_FCC") j=n*6;
		aux.Arw = new Histogram[j];
		aux.ArwC = new Histogram[j];
		for (i=0;i<j;i++) {
			aux.Arw[i].init(ether.hist1,ether.hist2,ether.histSteps);
			aux.ArwC[i].init(ether.hist1,ether.hist2,ether.histSteps);
		}
	}
	if (ether.isSet("optical")) aux.Sigma.init(0.0,2*ether.hist2,ether.histSteps);
	
	if (ether.isSet("ldos")) {
                aux.Ldos = new Histogram[n];
                for (i=0;i<n;i++) aux.Ldos[i].init(ether.hist1,ether.hist2,ether.histSteps);
        }
	
	aux.avMoments.insert(aux.avMoments.begin(),ether.tpem_cutoff,0.0);
	aux.curMoments.insert(aux.curMoments.begin(),ether.tpem_cutoff,0.0);
	
	aux.offdMoments.insert(aux.offdMoments.begin(),ether.tpem_cutoff*n,MatType(0.0,0.0));
	
	
	// Read files for dynvars if necessary
	string tmpString;
	vector<double> tmpVector;
	int color = int(ether.mpiRank/ether.mpiNop1);
	
	tmpVector.insert(tmpVector.begin(),n,0.0);
	
	if (ether.startType==4 || ether.startType==6) {
		if (ether.isSet("customconfig")) {
		 	std::stringstream ss;
        		std::string str;
        		ss << color;
        		ss >> str;
			ether.dynVarsInputFile=ether.dynVarsInputFile + str + ".sav";
		}
		tmpString=string("#Theta"); // + Ttoa(j);
		loadVector(tmpVector,ether.dynVarsInputFile,tmpString,ether.startLevel);
		if (ether.mpiRank==0) {
			cerr<<"Reading Theat from file "<<ether.dynVarsInputFile<<" at level "<<ether.startLevel;
			cerr<<" Theta[0]=="<<tmpVector[0]<<endl;
		}
		
		for (i=0;i<n;i++) {
			dynVars.theta[i]=tmpVector[i];
		}
		
		tmpString=string("#Phi"); // + Ttoa(j);
		loadVector(tmpVector,ether.dynVarsInputFile,tmpString,ether.startLevel);
		for (i=0;i<n;i++) {
			dynVars.phi[i]=tmpVector[i];
		}
#ifdef MODEL_KONDO_INF_TWOBANDS
		if (!ether.isSet("freezephonons")) {
			for (j=0;j<ether.D;j++) {
				tmpString=string("#Phonons")+ttos(j); // + Ttoa(j);
				loadVector(tmpVector,ether.dynVarsInputFile,tmpString,ether.startLevel);
				for (i=0;i<n;i++) dynVars.phonons[i][j]=tmpVector[i];
			}
					
		}
#endif
	}
	if (ether.startType==5) {
		int L = geometry.length();
		int column, row;
		for (row = 0; row < L; row++) {
			for (column=0;column<L;column++) {
				dynVars.theta[row + column*L]=M_PI*0.5;
			}
			if (row % 4 == 0)
				for (column=3; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
			
		
			if (row % 4 == 1) {
				for (column=1; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
				for (column=2; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
				for (column=3; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
			}
			
			if (row % 4 == 2)
				for (column=1; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
			
			if (row % 4 == 3) {
				for (column=0; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
				for (column=1; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
				for (column=3; column < L; column += 4)
					dynVars.phi[row*L+column]=M_PI;
				
			}
		}
		// if D==3 copy this plane onto the others:
		j= -1;
		if (geometry.dim()==3) {
			row=0;
			for (i=L*L;i<L*L*L;i++) {
				dynVars.theta[i]=dynVars.theta[row];
				dynVars.phi[i]=dynVars.phi[row];
				if (j!=1) dynVars.phi[i] += M_PI;
				row++;
				if (row==L*L) {
					row=0;
					j= -j;
				}
			}
		}
	}
}

// TPEM_FIXME: generalize
double calcNumber (vector<double> const &eig,Parameters const &ether,Aux &aux)
{
	unsigned int i;
	double beta=ether.beta,mu=aux.varMu,n_electrons;
	
	n_electrons=0;
	for (i=0;i<eig.size();i++) {
		n_electrons += utils::fermi((eig[i]-mu)*beta);
	}
	return n_electrons;
}


// TPEM_FIXME: generalize
double calcElectronicEnergy(vector<double>  const &eig,Parameters const &ether,Aux &aux)
{
	unsigned int i;
	double mu,beta,ee;
	
	mu =aux.varMu;
	beta = ether.beta;
	ee=0.0;
	
	for (i=0;i<eig.size();i++) {
		ee += eig[i]*utils::fermi((eig[i]-mu)*beta);
	}
	return ee;
		
}

double calcAction(vector<double>  const &eig,Parameters const &ether,Aux &aux)
{
	unsigned int i;
	double mu,beta,ee,tmp;
	
	mu =aux.varMu;
	beta = ether.beta;
	ee=0.0;
	
	for (i=0;i<eig.size();i++) {
		tmp=beta*(eig[i]-mu);
		ee += tmp+utils::logfermi(tmp);
		//cerr<<" calcaction "<<i<<" "<<ee<<" "<<tmp<<" "<<logfermi(tmp)<<endl;
	}
	return -ee;
		
}

void matrixPrint(MyMatrix<MatType> &matrix)
{
	int i,j,n=matrix.getRank();
	
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			cout<<i<<" "<<j<<" "<<real(matrix(j,i))<<" "<<imag(matrix(j,i))<<endl;
		}
	}
}


void setupHamiltonian(MyMatrix<MatType> &matrix,Geometry const &geometry,DynVars const &dynVars,
		Parameters const &ether,Aux &aux,int type)
{
	/*static int firstcall=1;
	static tpem_sparse *spMatrix;
	int matsize=ether.hilbertSize;
	
	if (firstcall) {
		spMatrix = new_tpem_sparse (matsize, ether.nonzero);
		firstcall=0;
	}
	
	kTpemHamiltonian (geometry,dynVars,spMatrix,ether,aux,type);*/
	createHamiltonian(geometry,dynVars,matrix,ether,aux,type);
	
	// no need to adjust spectrum since we're doing diagonalization here
	//crsToFull(matrix,spMatrix,ether,aux);

}

bool doMetropolis(vector<double> const &eNew,vector<double> const &eOld,
	Parameters const &ether, Aux &aux,double dsDirect,double sineupdate)
{
	unsigned int i;
	double X,beta,r,temp;
	double mu=aux.varMu;
	
	beta = ether.beta;
	X =1.0;
	
	for (i=0;i<eNew.size();i++) {
		
         if (eNew[i]>mu)
         	temp = (double)(1.0+exp(-beta*(eNew[i]-mu)))/(1.0+exp(-beta*(eOld[i]-mu)));
         else
             temp =(double)(1.0+exp(beta*(eNew[i]-mu)))/
					 (exp(beta*(eNew[i]-mu))+exp(-beta*(eOld[i]-eNew[i])));
         
		 X *= temp;
		 //cerr<<"temp="<<temp<<" X="<<X<<" "<<eNew[i]<<" "<<eOld[i];
		 //cerr<<" mu= "<<mu<<" i= "<<i<<endl;
     }
	
	if (ether.isSet("sineupdate")) X *= sineupdate;
	X *= (double) exp(-beta*dsDirect);
	X = (double)X/(1.0+X);
	
	r=ether.rng.myRandom(); 
	
//	cerr<<"X= "<<X<<" r="<<r<<endl;
	if (X>r) return true;
	else return false;
	
}


void r_newSpins(double thetaOld, double phiOld, double &thetaNew,double &phiNew, Parameters const &ether)
{
	if (ether.isSet("isingspins")) {
		if (thetaOld==0) thetaNew=M_PI; else thetaNew=0;
		phiNew=0;
	} else {
		if (ether.window<0) {
			thetaNew = 2*ether.rng.myRandom()-1;
			phiNew = 2*M_PI*ether.rng.myRandom();
			thetaNew = acos(thetaNew);
		} else {
			thetaNew=2*ether.rng.myRandom()- 1;
			if (thetaNew < -1) thetaNew= 0;
			if (thetaNew > 1) thetaNew = 0;		
			phiNew=phiOld+2*M_PI*(ether.rng.myRandom()- 0.5)*ether.window;
			thetaNew = acos(thetaNew);
		}
		if (ether.isSet("sineupdate")) {
			thetaNew = M_PI*ether.rng.myRandom();
		}
	
		
		
		while (thetaNew<0) {
			thetaNew = -thetaNew;
			phiNew+=M_PI;
		}	
		while (thetaNew>M_PI) {
			thetaNew -= M_PI;
			phiNew+=M_PI;
		}
			
		while (phiNew<0) phiNew += 2*M_PI;
		while (phiNew>2*M_PI) phiNew -= 2*M_PI;
	}
}

void r_newPhonons(vector<double> const &phononsOld,vector<double>  &phononsNew,Parameters const &ether)
{
	unsigned int i;
	for (i=0;i<phononsNew.size();i++) {
		phononsNew[i]=phononsOld[i] + (ether.rng.myRandom()- 0.5)*ether.window;
		//if (fabs(phononsNew[i]) > ether.maxPhonons) phononsNew[i]= 0.9*ether.maxPhonons;
	}
	
}

		
// TPEM_FIXME: moment manipulation instead of diagonalization
void doMonteCarlo(Geometry const &geometry,DynVars &dynVars,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	double phiNew,thetaNew,dsDirect=0,oldmu;
	bool flag;
	DynVars dynVars2;
	int dof;
	double sineupdate;
	vector<double> bcsPhiNew(2);
	
	size_t n = geometry.volume();
	
	vector<int> nac(ether.classFieldList.size(),0);
	vector<double> eigNewAllBands(ether.hilbertSize,0.0);
	dynVars2.theta=dynVars.theta;
	dynVars2.phi=dynVars.phi;
	dynVars2.phonons=dynVars.phonons;
	
	vector<double> phononsNew(ether.D,0.0);
	Phonons<Parameters,Geometry> phonons(ether,geometry);
	
	if (!ether.tpem) diag(aux.eigAllBands,geometry,dynVars2,ether,aux);
		
	if (ether.tpem && ether.carriers>0) calcMoments(dynVars,geometry,ether,aux,tpemOptions);
	
	
	for (size_t k=0; k<ether.classFieldList.size();k++) nac[k] = 0;
	
	for (size_t i=0;i<n;i++) {
		if (ether.modulus[i]==0) continue;
		
		for (size_t k=0;k<ether.classFieldList.size();k++) {
			dof = ether.classFieldList[k];
		
			if (ether.isSet("freezephonons") && dof==1) break;
			if (ether.isSet("freezespins") && dof==0) continue;
			
			if (ether.tpem) kTpemSetAllBands(geometry,dynVars,ether,aux);
		
			dynVars2.theta=dynVars.theta;
			dynVars2.phi=dynVars.phi;
			dynVars2.phonons=dynVars.phonons;
			
			switch (dof) {	
			case 0:    // spin dof
				r_newSpins(dynVars.theta[i],dynVars.phi[i],thetaNew,phiNew,ether);
				dynVars2.theta[i]=thetaNew;
				dynVars2.phi[i]=phiNew;		
				dsDirect = dSDirect(dynVars,dynVars2,i,geometry,ether);
				if (ether.isSet("tprime")) 
					dsDirect += directExchange2(dynVars2,geometry,ether)-directExchange2(dynVars,geometry,ether);
				break;
			case 1:  // phonons
				r_newPhonons(dynVars.phonons[i],phononsNew,ether);
				dynVars2.phonons[i]=phononsNew;
				dsDirect= dPhonons(dynVars,dynVars2,i,geometry,ether,phonons);
				break;
			}
			oldmu=aux.varMu;	
			if (ether.tpem) {
				aux.atmp=aux.varTpem_a;
				aux.btmp=aux.varTpem_b;
				if (ether.carriers>0) adjChemPotTpem(ether,aux,tpemOptions);
				flag=kTpemAllBands(i,geometry,dynVars2,ether,aux,dsDirect,tpemOptions);
			} else {
				diag(eigNewAllBands,geometry,dynVars2,ether,aux);
				if (ether.carriers>0) adjChemPot(eigNewAllBands,ether,aux);	
				sineupdate= sin(dynVars.theta[i]);
				if (dynVars.theta[i]!=0) {
					sineupdate = sin(dynVars2.theta[i])/sineupdate;
				} else {
					sineupdate = 1.0;
				}
				flag=doMetropolis(eigNewAllBands,aux.eigAllBands,ether,aux,dsDirect,sineupdate);					
			}
		
			if (flag && (ether.mcflag & 1)) { // Accepted
				switch (dof) {
				case 0:  // spin dof
					dynVars.theta[i]=thetaNew;
					dynVars.phi[i]=phiNew;
					break;
				case 1:
					dynVars.phonons[i]=phononsNew;
					break;
				}
				if (ether.tpem==0) copy(eigNewAllBands.begin(),eigNewAllBands.end(),aux.eigAllBands.begin());
				else tpem_sparse_copy(aux.sparseTmp[0],aux.sparseMatrix[0]);
																	
				nac[k]++;
			} else { // not accepted
				//aux.varMu=oldmu;
				if (ether.tpem && ether.isSet("adjusttpembounds")) { aux.varTpem_a=aux.atmp; aux.varTpem_b=aux.btmp; }
			}
		} // dof
				
	} // lattice sweep
	for (size_t k=0; k<ether.classFieldList.size();k++) aux.nac[k] += nac[k];

}



// TPEM_FIXME: ok, no changes needed
double calcMag(DynVars const &dynVars,Parameters const &ether)
{
	int i;
	int n =ether.linSize;
	vector<double> mag(3);
	
	mag[0]=mag[1]=mag[2]=0.0;
	for (i=0;i<n;i++) {
		if (ether.modulus[i]==0) continue;
		mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
		mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
		mag[2] += cos(dynVars.theta[i]);
	}
	return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
	
}

// TPEM_FIXME: ok, no changes needed
void calcClasCor(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,vector<double> &cc,
        vector<double> &weight)
{
	int i,j,k;
	int n = geometry.volume();
	double temp;
	int counter;
	
	for (i=0;i<n;i++) {
		temp=0.0;
		counter=0;
		for (j=0;j<n;j++) {
			k = geometry.add(i,j);
			//cerr<<"calcClasCor: "<<i<<"+"<<j<<"="<<k<<endl;
			if (ether.modulus[k]==0 || ether.modulus[j]==0) continue;
			temp+= cos(dynVars.theta[k])*cos(dynVars.theta[j])+
					sin(dynVars.theta[k])*sin(dynVars.theta[j])*
					cos(dynVars.phi[k]-dynVars.phi[j]);
			counter++;
		}
		if (counter>0) temp /= counter;
		weight[i]=counter;
		cc[i] = temp;
	}
}

double calcClasCor(DynVars const &dynVars,Parameters const &ether, int i,int j)
{
	double temp=0;
	if (ether.modulus[i]==0 || ether.modulus[j]==0) return 0;
	temp= cos(dynVars.theta[i])*cos(dynVars.theta[j])+
		sin(dynVars.theta[i])*sin(dynVars.theta[j])*cos(dynVars.phi[i]-dynVars.phi[j]);
	return temp;
}


double calcSq(vector<double> const &correlation,Geometry const &geometry,int q)
{
        int i,l=geometry.length(0),n=geometry.volume();
        double scprod,sum=0;

        for (i=0;i<n;i++) {
                scprod = geometry.scalarProd(i,q);
                scprod *= (double)(2.*M_PI/l);
                sum += correlation[i] * scprod;
        }
        return sum;

}

void rotateSpins(DynVars &dynVars1,DynVars const &dynVars2, Geometry const &geometry,int dir)
{
	if (!geometry.isCubicType() && !(geometry.latticeName()=="prism") && !(geometry.latticeName()=="rectangle")) {
		cerr<<"rotateSpins: not supported for lattice type="<<geometry.latticeName()<<endl;
		cerr<<"Exiting at this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1); //ERROR CODE		
	}
	
	int lx,ly,lz,x,y,z,d = geometry.dim();
	
	switch (d) {
	case 1:
		cerr<<"Rotation of a one-dimensional lattice is strictly forbidden\n";
		cerr<<"Exiting at this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1); //ERROR CODE		
		break;
	case 2:
		lx=geometry.length(0); ly=geometry.length(1); lz=geometry.length(2);
		for (y=0;y<ly;y++) {
			for (x=0;x<lx;x++) {
				copyDynVars(dynVars1 ,x+y*lx,dynVars2,ly-1-y+lx*x);
			}
		}
		break;
	case 3:
		lx=geometry.length(0); ly=geometry.length(1); lz=geometry.length(2);
		if (dir==0) { // FIXME: ADD THE THREE DIRECTIONS
			for (z=0;z<lz;z++) {
				for (y=0;y<ly;y++) {
					for (x=0;x<lx;x++) {
						copyDynVars(dynVars1 ,x+y*lx+z*lx*ly,dynVars2,x+lx*(lz-1-z)+y*lx*ly);
					}
				}
			}
		} else {
			for (z=0;z<lz;z++) {
				for (y=0;y<ly;y++) {
					for (x=0;x<lx;x++) {
						copyDynVars(dynVars1 ,x+y*lx+z*lx*ly,dynVars2,lz-1-z+lx*y+x*lx*ly);
					}
				}
			}
		}
		break;
	}
}

double Density(double x)
{
	double	ret;
	double	a, b, mu, beta;
	
	tmpValues(a,b,mu,beta,1);
	ret = 0.5 * (1.0 - tanh (0.5 * beta * (a * x + b - mu)));
	return ret; 
}

double effective_action (double x)
{
	double  a, b, mu, beta , ret;
	
	
	tmpValues(a,b,mu,beta,1);
	
	if (a * x + b >= mu)
		ret =	log (1.0 + exp (-beta * (a * x + b - mu)));
	else
		ret =	-beta * (a * x + b - mu)
			+ log (1.0 + exp (beta * (a * x + b - mu)));
	
	return   ret;
}

double Energy (double x)
{
	double	ret;
	double	a, b, mu, beta;
	
	
	tmpValues(a,b,mu,beta,1);
	ret = (a * x + b) * 0.5 * (1.0 - tanh (0.5 * beta * (a * x + b - mu)));
	
	return ret;
}

double measure_density  (vector<double> const &moment,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff = ether.tpem_cutoff;
	static vector<double>	coeffs(cutoff);
	double ret;
	static int firstcall=1;
	double beta=ether.beta;

	if (firstcall || ether.carriers>0) {
		firstcall=0;
		tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
		tpem_calculate_coeffs (coeffs, Density,tpemOptions);	
	}
	tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
	if (ether.isSet("adjusttpembounds")) tpem_calculate_coeffs (coeffs, Density,tpemOptions);		
	
	ret = tpem_expansion (moment, coeffs);
	return ret;
}

double measure_energy (vector<double> const &moment,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff = ether.tpem_cutoff;
	static vector<double>	coeffs(cutoff);
	double ret;
	static int firstcall=1;
	double beta=ether.beta;
	
	if (firstcall || ether.carriers>0) {
		firstcall=0;
		tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
		//cout<<"measure_energy: Calculating Coefficients\n";
		tpem_calculate_coeffs (coeffs, Energy,tpemOptions);		
	}
	tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
	if (ether.isSet("adjusttpembounds")) tpem_calculate_coeffs (coeffs, Energy,tpemOptions);
		
	
	ret = tpem_expansion (moment, coeffs);
	return ret;
}

double measure_action (vector<double> const &moment,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	
	double ret;
	int cutoff = ether.tpem_cutoff;
	static vector<double>	coeffs(cutoff);
	static int firstcall=1;
	double beta=ether.beta;
	
	
	if (firstcall || ether.carriers>0) {
		firstcall=0;
		tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
		//cout<<"measure_action: Calculating Coefficients\n";
		tpem_calculate_coeffs (coeffs, effective_action,tpemOptions);
	}
	tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
	if (ether.isSet("adjusttpembounds")) tpem_calculate_coeffs (coeffs, effective_action,tpemOptions);
		
	ret = tpem_expansion (moment, coeffs);
	return ret;
}



void calcMoments(vector<double> &moment,DynVars const &dynVars,Geometry const &geometry,Parameters const &ether,
	Aux &aux,TpemOptions const &tpemOptions)
{	
	
	kTpemHamiltonian (geometry, dynVars,aux.sparseMatrix[0],ether,aux,0);
	if (ether.isSet("adjusttpembounds")) {
		if (tpemAdjustBounds(aux.sparseMatrix[0],ether,aux)!=0) {
			cerr<<"Cannot adjust bounds for tpem spectrum\n";
			exit(1);
		}
	}
	tpem_calculate_moment (aux.sparseMatrix[0], moment, tpemOptions);
	//if (ether.mpiRank==0) cerr<<"moment[0]="<<moment[0]<<endl;
}

void calcMoments(DynVars const &dynVars,Geometry const &geometry,Parameters const &ether,Aux &aux,TpemOptions const
&tpemOptions)
{
	vector<double> moment(ether.tpem_cutoff);
	calcMoments(moment,dynVars,geometry,ether,aux,tpemOptions);
	for (unsigned int i=0;i<moment.size();i++) aux.curMoments[i] = moment[i];
	
}



			
// TPEM_FIXME: Electronic observables must be modified
template<typename ConcurrencyIoType>
void doMeasurements(int iter,DynVars const &dynVars,Geometry const &geometry,Io<ConcurrencyIoType> &io,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	double n_electrons;
	unsigned int n=geometry.volume();
	double temp,temp2;
	int cutoff = ether.tpem_cutoff;
	vector<double> moment(cutoff);
	string s;
	int d = geometry.dim();
	DynVars dynVars2;
	
	vector<double> tmp(n,0),tmpweight(n,0);
	
	s = "iter=";
	io.historyPrint(s,iter);
	
	n_electrons=0;
		
	if (ether.tpem) {
		calcMoments(moment,dynVars,geometry,ether,aux,tpemOptions);
		n_electrons=measure_density (moment, ether,aux,tpemOptions);
	} else {
	
		diag(aux.eigOneBand,geometry,dynVars,ether,aux);
		
		for (int i=0;i<ether.hilbertSize;i++) {
			//aux.eigM[i] += aux.eigOneBand[i];
		}
#ifdef MODEL_KONDO_DMS_MANYBANDS
		for (int i=0;i<ether.numberOfOrbitals;i++) {
			temp = calcNumber(dynVars,geometry,ether,aux,i); // number of electrons for band i  
			s ="Number_Of_Band"+ttos(i)+"=";
			io.historyPrint(s,temp);
		}
#endif
		n_electrons=calcNumber (aux.eigOneBand,ether,aux);

	}
		
		
	s ="Number_Of_Electrons=";
	io.historyPrint(s,n_electrons);
	
	s = "rankGlobal=";
	double varRank = ether.mpiRank;
	io.historyPrint(s,varRank);
	
		
	if (ether.tpem) temp=measure_energy (moment, ether,aux,tpemOptions);
	else temp=calcElectronicEnergy(aux.eigOneBand,ether,aux);
	
	s="Electronic Energy=";
	io.historyPrint(s,temp);
	temp2=calcSuperExchange(dynVars,geometry,ether);
	s="Superexchange=";
	io.historyPrint(s,temp2);
	temp += temp2;
	if (ether.isSet("tprime")) {
		temp2=directExchange2(dynVars,geometry,ether);
		s="Superexchange2=";
		io.historyPrint(s,temp2);
		temp += temp2;
	}
#ifdef MODEL_KONDO_INF_TWOBANDS	
	Phonons<Parameters,Geometry> phonons(ether,geometry);	
	temp2 = calcPhononEnergy(dynVars,geometry,ether,phonons);
	s="PhononEnergy=";
	io.historyPrint(s,temp2);
	temp += temp2;
#endif		
	// total energy = electronic energy + superexchange + phonon energy
	s="TotalEnergy=";
	io.historyPrint(s,temp);
		
	if (ether.tpem) temp=measure_action (moment, ether,aux,tpemOptions);
	else temp=calcAction(aux.eigOneBand,ether,aux);
	s="Action=";
	io.historyPrint(s,temp);
	
	s="Number_Of_Holes=";
	io.historyPrint(s,ether.hilbertSize-n_electrons);

	s="ChemPot=";
	io.historyPrint(s,aux.varMu);
	
	if (ether.isSet("adjusttpembounds") || iter==0) {
		s="tpem_a=";
		io.historyPrint(s,aux.varTpem_a);
		s="tpem_b=";
		io.historyPrint(s,aux.varTpem_b);
	}
	
	for (size_t i=0;i<n;i++) {
		tmp[i]=0.0;
		tmpweight[i]=0.0;
	}
	calcClasCor(geometry,dynVars,ether,tmp,tmpweight);
	if (iter==0 && ether.mpiRank==0 && ether.conc<ether.linSize) {
		cout<<"# weight = ";
		for (size_t i=0;i<n;i++) {
			cout<<tmpweight[i]<<" ";
		}
		cout<<endl;
	}
	
	s="Cor[1]=";
	io.historyPrint(s,tmp[1]);

	if (geometry.isCubicType()) {
		int l = geometry.length();
		int i=0;
		switch (d) {
			case 1:
				i=int(l*0.5);
			case 2:
				i+=int(l*l*0.5);
			case 3:
				i+=int(l*l*l*0.5);
		}
		temp=calcSq(tmp,geometry,i);
       
		s="Sq[pi]=";
		io.historyPrint(s,temp);
	
		if (d==2) {
			i=l*0.5;
			temp=calcSq(tmp,geometry,i);
			s="Sq[(pi,0)]=";
			io.historyPrint(s,temp);
        	}
	}
	for (size_t i=0;i<n;i++) aux.clasCor[i] += tmp[i];
	
	if (iter==0) {
		s="Weight[0]=";
		io.historyPrint(s,tmpweight[0]);
	}
	temp= calcClasCor(dynVars,ether,0,4);
	s="CustomCorrelation=";
        io.historyPrint(s,temp);
	
	temp = calcMag(dynVars,ether);
	s="Mag2=";
	io.historyPrint(s,temp);
	
	kTpemMoments(moment,aux,ether);
	
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
	temp=calcPhononV(geometry,dynVars,ether);
	temp2=calcPhononV2(geometry,dynVars,ether);
	s="calcPhononV=";
	io.historyPrint(s,temp);
	s="calcPhononV2=";
	io.historyPrint(s,temp2);
	temp=temp2/ether.linSize - square(temp/ether.linSize);
	s="phononDeltaV=";
	io.historyPrint(s,temp);	
#endif		
	for (size_t dof=0;dof<ether.classFieldList.size();dof++) {
		aux.nac[dof] = (double)aux.nac[dof]/ether.iterUnmeasured;
		s ="Accepted["+ttos(dof)+"]=";
		io.historyPrint(s,aux.nac[dof]);
		aux.nac[dof]=0;
	}
	
	if (!ether.tpem) {
		diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
		
		accNOfOmega(geometry,dynVars,ether,aux);
		accLcd(geometry,dynVars,ether,aux);
		if (ether.isSet("savelcd")) io.printLcd(ether,aux);
		temp=calcKinetic(dynVars,geometry,ether,aux);
		if (ether.isSet("akw")) accAkw(geometry,dynVars,ether,aux);
		if (ether.isSet("optical")) accOptical(geometry,dynVars,ether,aux);
		if (ether.isSet("ldos")) accLdos(geometry,dynVars,ether,aux);	
		if (ether.isSet("chargecorrelation")) accChargeCorrelation(0,0,geometry,dynVars,ether,aux);
		if (ether.isSet("chargecorrelation")) accChargeCorrelation(0,1,geometry,dynVars,ether,aux);
		if (ether.isSet("chargecorrelation")) accChargeCorrelation(1,0,geometry,dynVars,ether,aux);
		if (ether.isSet("chargecorrelation")) accChargeCorrelation(1,1,geometry,dynVars,ether,aux);
		if (ether.isSet("orbitalangles")) accOrbitalAngles(geometry,ether,aux);
		if (ether.isSet("orbitalcorrelation")) accOrbitalCorrelation(0,0,geometry,dynVars,ether,aux);
		if (ether.isSet("orbitalcorrelation")) accOrbitalCorrelation(0,1,geometry,dynVars,ether,aux);
		if (ether.isSet("orbitalcorrelation")) accOrbitalCorrelation(1,0,geometry,dynVars,ether,aux);
		if (ether.isSet("orbitalcorrelation")) accOrbitalCorrelation(1,1,geometry,dynVars,ether,aux);
		//std::vector<size_t> q(3);
		//q[0]=0;
		//q[1]=128;
		//q[2]=129;
		if (ether.isSet("nanocluster")) {
			psimag::Matrix<std::complex<double> > sq(ether.q.size(),ether.linSize);
	//				ether.linSize);
			calcLocalk(sq,ether.q,geometry,dynVars,ether);
			io.historyPrint("nanocluster",sq);
		}
	
	} else {
		temp=measure_kinetic(geometry,ether,aux,tpemOptions);
		if (ether.isSet("akw")) kTpemMomentsCl(geometry,ether,aux,tpemOptions);
		if (ether.isSet("optical")) kTpemMomentsOptical(geometry,ether,aux,tpemOptions);
		
	}


	s ="KineticEnergy=";
	io.historyPrint(s,temp);

	if (ether.isSet("conductance") &&   ether.mpiTpemRank==0  ) {
		int maxiter = 1000;
		double maxerror = 1.e-6;
		
		setupHamiltonian(aux.matrix,geometry,dynVars,ether,aux,0);
		temp=conductance3d(aux.matrix,aux.varMu,ether.hilbertSize,geometry.dim(),geometry.length(),maxiter,maxerror);
		copyDynVars(dynVars2,dynVars);
		// cerr<<"temp="<<temp<<endl; CONDUCTANCEDIRORIGINAL
		if (geometry.dim()>1) {
			rotateSpins(dynVars2,dynVars,geometry,0);
			setupHamiltonian(aux.matrix,geometry,dynVars2,ether,aux,0);
			temp2=conductance3d(aux.matrix,aux.varMu,ether.hilbertSize,geometry.dim(),geometry.length(),maxiter,maxerror);
			// cerr<<"temp2="<<temp2<<endl; CONDUCTANCEDIR0
			temp += temp2;
		}
		if (geometry.dim()>2) {
			rotateSpins(dynVars2,dynVars,geometry,1);
			setupHamiltonian(aux.matrix,geometry,dynVars2,ether,aux,0);
			temp2=conductance3d(aux.matrix,aux.varMu,ether.hilbertSize,geometry.dim(),geometry.length(),maxiter,maxerror);
			// cerr<<"temp2="<<temp2<<endl; CONDUCTANCEDIR1
			temp += temp2;
		}
		temp /= geometry.dim();
		s="Conductance=";
		io.historyPrint(s,temp);
	}
	
	if (ether.isSet("saveall")) io.printSnapshot(dynVars,ether);
	if (ether.isSet("filelastconfig")) io.printSnapshot(dynVars,ether,1);

}				
