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
#include "argument.h"
#include "conductance.h"
#include "VectorGenerator.h"
#include "RandomGenerator.h"
#include "MpiParameter.h"

bool Io::isInstatiated=false;
bool Io::isInit=false;

using namespace std;

extern void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type);
extern void setupHamiltonian(MyMatrix<MatType> & matrix,Geometry const &geometry, DynVars const &dynVars, 
        Parameters const &ether,Aux &aux,int bandindex);
extern void setHilbertParams(Parameters &ether, Aux &aux, Geometry const &geometry);
extern void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry);


void setTheRankVector(Parameters& ether,std::vector<size_t>& v,std::vector<size_t>& w)
{
	
	size_t size0 = ether.numberOfBetas;
	v.resize(2);
	w.resize(2);
	v[0] = ether.mpiRank % size0;
	v[1] = size_t(ether.mpiRank/size0);
	w[0] = size0;
	w[1] = ether.mpiSize / size0;
}

void registerHook(Parameters& ether)
{

	setTheRankVector(ether,ether.localRank,ether.localSize);
	//example of non-random
	VectorGenerator<double> betaGenerator(ether.betaVector,ether.localRank[0]);
	typedef MpiParameter<double,VectorGenerator<double>,Parameters > MpiParameterBeta;
	MpiParameterBeta beta(ether.beta,ether,betaGenerator,MpiParameterBeta::SEPARATE,ether.localRank[0]); // beta 4 10 20 30 40
	
	// example of random
	//typedef MpiParameter<std::vector<double>,RandomGenerator,Parameters> MpiParameterJaf;
	//RandomGenerator jafGenerator("bimodal",2,3,ether.localSize[1],ether.localRank[1]); // jaf first, deltaJAf second
	//MpiParameterJaf jafvector(ether.JafVector,ether,jafGenerator,MpiParameterJaf::TOGETHER,ether.localRank[1]);

}

int spf_entry(int argc,char *argv[])
{
	Parameters ether;
	char sfile[256];
	TpemOptions tpemOptions;
	
	tpemOptions.mpi_nop2 =1;
	
	rlimit myrlim;
        myrlim.rlim_cur=myrlim.rlim_max=0;
        setrlimit(RLIMIT_CORE, &myrlim);

        tpem_init(argc,argv,tpemOptions,sfile);
	ether.mpiRank = tpemOptions.rank;
	ether.mpiNop1 = tpemOptions.mpi_nop1;
	ether.mpiNop2 = tpemOptions.mpi_nop2;
	ether.mpiTpemRank = tpemOptions.tpem_rank;
	ether.mpiSize = ether.mpiNop1  * ether.mpiNop2;
	
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
	Io io;
	
	srand(time(0));

	
	if (io.input(sfile,geometry,dynVars,ether,aux)) {
		if (ether.mpiRank==0) cerr<<"A problem ocurred while trying to load file "<<sfile<<endl;
		return 1;
	}

	registerHook(ether);	
	
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
		printProgress(iter,ether.iterTherm,10,'*',ether.mpiRank);
		doMonteCarlo(geometry,dynVars,ether,aux,tpemOptions);
	}
	if (ether.mpiRank==0) {
		io.fileOut()<<"\nMEASUREMENT PHASE ----------\n";
		io.fileOut()<<"PROGRESS:            "<<endl;
		cout.flush();
	}

	for (iter=0;iter<ether.iterEffective;iter++) {
		printProgress(iter,ether.iterEffective,10,'*',ether.mpiRank);

		
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
                ether.jafFile = ether.jafFile + dtos(y) + ".txt";
                string s2="#Jafvector";
                loadVector(ether.JafVector,ether.jafFile,s2,1);
                ether.options = ether.options + ",jafvector";
        }
	
}

double calcBcsDelta2(DynVars const &dynVars,Parameters const &ether)
{
	unsigned int i;
	double s=0;
	for (i=0;i<dynVars.bcsDelta.size();i++) {
		s += square(dynVars.bcsDelta[i]);
	}
	return s;
}

double calcBcsDeltaV2(DynVars const &dynVars,Parameters const &ether)
{
	unsigned int i,j;
	double s=0;
	double tmp=0;
	
	for (i=0;i<dynVars.bcsDelta.size();i++) {
		tmp=0;
		for (j=0;j<ether.D;j++) tmp += ether.bcsV[i+j*ether.linSize];
		tmp /= ether.D;
		s += tmp * square(dynVars.bcsDelta[i]);
	}
	return s;
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
        double tmp;

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



double calcPhononDiff(int direction,int ind,DynVars const &dynVars,Geometry const &geometry)
{
	if (direction >= geometry.dim()) return 0; 
	int j = geometry.neighbor(ind,2*direction+1);
	return  (dynVars.phonons[ind][direction]-dynVars.phonons[j][direction]);
}
		
double calcPhonon(int ind,DynVars const &dynVars,Geometry const &g,Parameters const &ether,int what)
{
	double ret=0;
	double sqrt3=1.732050807569;
	double sqrt2=1.414213562373;
	double sqrt6=2.449489742783;
	
	if (what==0)  { /* calc q1 */
		ret = (calcPhononDiff(0,ind,dynVars,g) + calcPhononDiff(1,ind,dynVars,g) +
		 calcPhononDiff(2,ind,dynVars,g));
		ret /= sqrt3;
	} else if (what==1) { /* calc q2 */
		ret = (calcPhononDiff(0,ind,dynVars,g)-calcPhononDiff(1,ind,dynVars,g));
		ret /= sqrt2;
	} else if (what==2) { /* calc q3 */
		ret = (2*calcPhononDiff(2,ind,dynVars,g)-calcPhononDiff(0,ind,dynVars,g)-calcPhononDiff(1,ind,dynVars,g));
		ret /= sqrt6;
	} else {
		cerr<<"I don't know what to calculate at "<<__FILE__<<" "<<__LINE__<<" with what="<<what<<endl;
		exit(1);
	}
	return ret;
}

double dPhonons(DynVars const &dynVars,DynVars const &dynVars2,int i,Geometry const &geometry, Parameters const &ether)
{
	double tmp;
	int alpha,k,j;
	double dE=0.0;
	
#ifdef MODEL_KONDO_INF_TWOBANDS
	for (alpha=0;alpha<dynVars.phonons[i].size();alpha++) {
		tmp=0;
		tmp += square(calcPhonon(i,dynVars2,geometry,ether,alpha));
		tmp -= square(calcPhonon(i,dynVars,geometry,ether,alpha));
		for (k=0;k<geometry.z(i);k++) {
			j = geometry.neighbor(i,k);
			tmp += square(calcPhonon(j,dynVars2,geometry,ether,alpha));
			tmp -= square(calcPhonon(j,dynVars,geometry,ether,alpha));
		}
		dE += ether.phononEd[alpha]*tmp;
	}
#else	
	for (j=0;j<dynVars.phonons[i].size();j++)
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
		dS += ether.jafprime*tmp;
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

double calcPhononEnergy(DynVars const &dynVars, Geometry const &geometry, Parameters const &ether)
  // Written by IS Oct-19-04.
  // Calculates pure phononic energy of the system (last term in the Hamiltonian).
{
  int n = geometry.volume();
  double tmp;
  int i,j;
  unsigned int alpha;
  double dE=0.0;
	
  for (i=0;i<n;i++) 
    {

#ifdef MODEL_KONDO_INF_TWOBANDS
      for (alpha=0;alpha<dynVars.phonons[i].size();alpha++) 
	{
	  tmp = calcPhonon(i,dynVars,geometry,ether,alpha);
	  // cerr << "site "<< i << "phonon "<< alpha << ": " <<tmp << endl;
	  dE += ether.phononEd[alpha]*square(tmp);
	}
#else	
      for (j=0;j<dynVars.phonons[i].size();j++)
	dE += square(dynVars.phonons[i][j]);
#endif
    }
  return dE*0.5;
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
	double tmp2,realpart,imagpart;
	
	n=a.getRank();
	
	for (i = 0; i < n ; i++) for (k=0;k<n;k++) a.set(i,k,0.0);
		
	for (i = 0; i < n; i++){
		for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++){
			//a[i * n + matrix->colind[k]] = matrix->values[k];
			
			
			if(i==matrix->colind[k]){
				tmp2=aux.varTpem_b;
			}
			else{
				tmp2=0;
			}
			realpart=real(matrix->values[k])*aux.varTpem_a+tmp2;
			imagpart=imag(matrix->values[k])*aux.varTpem_a;
			a.set(matrix->colind[k],i,MatType(realpart,imagpart));
		}
	}
	if (!a.isHermitian()) {
		std::cerr<<a<<endl;
		exit(1);
	}
	
}

double calculate_dS (vector<double> &moment,tpem_sparse *matrix1, tpem_sparse *matrix, 
	vector<unsigned int> const &support,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff = ether.tpem_cutoff;
        static vector<double>  coeffs(cutoff),moment0(cutoff),moment1(cutoff);
        double dS;
        static int firstcall=1;
        double a=aux.varTpem_a,b=aux.varTpem_b,mu=aux.varMu,beta=ether.beta;
        int ii;
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
	int i,j;
	if (eig.size() != ether.hilbertSize) {
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
		int blocksize=aux.matrix.getRank()/2;
		vector<double> eig1(blocksize),eig2(blocksize);
		
		
		lowerMatrix.init(blocksize,0);
		upperMatrix.init(blocksize,0);
		for (i=0;i<blocksize;i++) {
			for (j=0;j<blocksize;j++)  {
				upperMatrix(i,j)=aux.matrix(i,j);
				lowerMatrix(i,j)=aux.matrix(i+blocksize,j+blocksize);
			}
		}
		diag(upperMatrix,eig1,jobz);
		diag(lowerMatrix,eig2,jobz);
		for (i=0;i<blocksize;i++) {
			eig[i]=eig1[i];
			eig[i+blocksize]=eig2[i];
		}
		if (jobz=='V') { // eigenvectors
			for (i=0;i<blocksize;i++) {
				for (j=0;j<blocksize;j++) {
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
		for (i=0;i<eig.size();i++) {
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
	unsigned int j;
	

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
		if (dS>myRandom()) accept=true;
		else accept=false;
	}	else {
		accept = (dS > 0.0 || myRandom () < exp (dS));
	}
	if (accept) {
		// update current moments
		for (j=0;j<ether.tpem_cutoff;j++) aux.curMoments[j] -= moment[j];		
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
	
	if (ether.isSet("randomize")) myRandomSeed(time(0));
	
	// NOTE: THIS FUNCTION SHOULD NOT BE REENTRANT--> FIXME
		
	setHilbertParams(ether,aux,geometry);
	
	aux.nac.insert(aux.nac.begin(),ether.classFieldList.size(),0);
	
	dynVars.theta.insert(dynVars.theta.begin(),n,0.0);
	dynVars.phi.insert(dynVars.phi.begin(),n,0.0);
	
	dynVars.Tx.insert(dynVars.Tx.begin(),n,0.0);
	dynVars.Ty.insert(dynVars.Ty.begin(),n,0.0);
	dynVars.Tz.insert(dynVars.Tz.begin(),n,0.0);
	
	for (i=0;i<d;i++) vtmp.push_back(0.0);
	
	if (ether.bcsDelta0>0) {
		dynVars.bcsDelta.insert(dynVars.bcsDelta.begin(),n,ether.bcsDelta0);
		for (i=0;i<n;i++) dynVars.bcsPhi.push_back(vtmp);
	}
	
	
	for (i=0;i<n;i++) {
		dynVars.phonons.push_back(vtmp);
		
	}
	for (i=0;i<n;i++) {
		switch (ether.startType) {
			case 1:
				dynVars.theta[i]=M_PI*myRandom();
				dynVars.phi[i]=2*M_PI*myRandom();
				if (ether.bcsDelta0>0) {
					//dynVars.bcsDelta[i] = ether.bcsDelta0*myRandom();
					for (j=0;j<dynVars.bcsPhi[i].size();j++) dynVars.bcsPhi[i][j] = 2*M_PI*myRandom();
				}
				if (!ether.isSet("freezephonons")) {
					if (ether.isSet("verbose") && ether.mpiRank==0) cerr<<"Randomizing phonons\n";
					for (j=0;j<d;j++) dynVars.phonons[i][j]=ether.maxPhonons*(myRandom()-0.5);
				}								
				break;
			case 2:
				dynVars.theta[i]=(parity(i,d,geometry.length())==1) ? 0 : M_PI;
				break;
			case 3:
				dynVars.theta[i]=M_PI*0.5;
				dynVars.phi[i]=M_PI*0.5;
				if (ether.bcsDelta0>0) {
					//for (j=0;j<dynVars.bcsPhi[i].size();j++) dynVars.bcsPhi[i][j] = M_PI*0.5;
					dynVars.bcsPhi[i][0] = 0;
					dynVars.bcsPhi[i][1] = M_PI;
				}
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
		n_electrons += fermi((eig[i]-mu)*beta);
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
		ee += eig[i]*fermi((eig[i]-mu)*beta);
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
		ee += tmp+logfermi(tmp);
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
	static int firstcall=1;
	static tpem_sparse *spMatrix;
	int matsize=ether.hilbertSize;
	
	if (firstcall) {
		spMatrix = new_tpem_sparse (matsize, ether.nonzero);
		firstcall=0;
	}
	
	kTpemHamiltonian (geometry,dynVars,spMatrix,ether,aux,type);
	
	// no need to adjust spectrum since we're doing diagonalization here
	crsToFull(matrix,spMatrix,ether,aux);

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
	
	r=myRandom(); 
	
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
			thetaNew = 2*myRandom()-1;
			phiNew = 2*M_PI*myRandom();
			thetaNew = acos(thetaNew);
		} else {
			thetaNew=2*myRandom()- 1;
			if (thetaNew < -1) thetaNew= 0;
			if (thetaNew > 1) thetaNew = 0;		
			phiNew=phiOld+2*M_PI*(myRandom()- 0.5)*ether.window;
			thetaNew = acos(thetaNew);
		}
		if (ether.isSet("sineupdate")) {
			thetaNew = M_PI*myRandom();
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
		phononsNew[i]=phononsOld[i] + (myRandom()- 0.5)*ether.window;
		//if (fabs(phononsNew[i]) > ether.maxPhonons) phononsNew[i]= 0.9*ether.maxPhonons;
	}
	
}

void r_newBcsFields(double bcsDeltaOld,vector<double> const &bcsPhiOld,double &bcsDeltaNew,vector<double> &bcsPhiNew,Parameters const
&ether)
{
	unsigned int i;
	bcsDeltaNew = bcsDeltaOld; //ether.bcsDelta0 * myRandom();
	for (i=0;i<bcsPhiNew.size();i++) bcsPhiNew[i] = bcsPhiOld[i]+ M_PI * (0.5 -myRandom());
}
		
// TPEM_FIXME: moment manipulation instead of diagonalization
void doMonteCarlo(Geometry const &geometry,DynVars &dynVars,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int i,n,k;
	double phiNew,thetaNew,dsDirect=0,oldmu;
	bool flag;
	DynVars dynVars2;
	int dof;
	double sineupdate,bcsDeltaNew;
	vector<double> bcsPhiNew(2);
	
	n = geometry.volume();
	
	vector<int> nac(ether.classFieldList.size(),0);
	vector<double> eigNewAllBands(ether.hilbertSize,0.0);
	dynVars2.theta=dynVars.theta;
	dynVars2.phi=dynVars.phi;
	dynVars2.phonons=dynVars.phonons;
	if (ether.bcsDelta0>0) {
		dynVars2.bcsDelta=dynVars.bcsDelta;
		dynVars2.bcsPhi=dynVars.bcsPhi;
	}
	
	vector<double> phononsNew(ether.D,0.0);
	
	
	if (!ether.tpem) diag(aux.eigAllBands,geometry,dynVars2,ether,aux);
		
	if (ether.tpem && ether.carriers>0) calcMoments(dynVars,geometry,ether,aux,tpemOptions);
	
	
	for (k=0; k<ether.classFieldList.size();k++) nac[k] = 0;
	
	for (i=0;i<n;i++) {
		if (ether.modulus[i]==0) continue;
		
		for (k=0;k<ether.classFieldList.size();k++) {
			dof = ether.classFieldList[k];
		
			if (ether.isSet("freezephonons") && dof==1) break;
			if (ether.isSet("freezespins") && dof==0) continue;
			
			if (ether.tpem) kTpemSetAllBands(geometry,dynVars,ether,aux);
		
			dynVars2.theta=dynVars.theta;
			dynVars2.phi=dynVars.phi;
			dynVars2.phonons=dynVars.phonons;
			if (ether.bcsDelta0>0) {
				dynVars2.bcsDelta=dynVars.bcsDelta;
				dynVars2.bcsPhi=dynVars.bcsPhi;
			}
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
				dsDirect= dPhonons(dynVars,dynVars2,i,geometry,ether);
				break;
			case 2: // bcs delta/phi
				r_newBcsFields(dynVars.bcsDelta[i],dynVars.bcsPhi[i],bcsDeltaNew,bcsPhiNew,ether);
				dynVars2.bcsDelta[i] = bcsDeltaNew;
				dynVars2.bcsPhi[i] = bcsPhiNew;
				dsDirect = (calcBcsDeltaV2(dynVars2,ether)-calcBcsDeltaV2(dynVars,ether)); 
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
				case 2:
					dynVars.bcsDelta[i] = bcsDeltaNew;
					dynVars.bcsPhi[i] = bcsPhiNew;
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
	for (k=0; k<ether.classFieldList.size();k++) aux.nac[k] += nac[k];

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

// TPEM_FIXME: ok, no changes needed
void calcBcsPhiCor(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,vector<double> &cc,
	int link1,int link2)
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
			temp += cos(dynVars.bcsPhi[j][link1]-dynVars.bcsPhi[k][link2]);
		}
		cc[i] += temp;
	}
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
void doMeasurements(int iter,DynVars const &dynVars,Geometry const &geometry,Io &io,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	double n_electrons;
	unsigned int i;
	unsigned int n=geometry.volume();
	double temp,temp2;
	int cutoff = ether.tpem_cutoff;
	vector<double> moment(cutoff);
	int dof;
	string s;
	int d = geometry.dim();
	DynVars dynVars2;
	
	vector<double> tmp(n,0),tmpweight(n,0);
	
	s = "iter="+ttos(iter);
	io.historyPrint(s);
	
	n_electrons=0;
		
	if (ether.tpem) {
		calcMoments(moment,dynVars,geometry,ether,aux,tpemOptions);
		n_electrons=measure_density (moment, ether,aux,tpemOptions);
	} else {
	
		diag(aux.eigOneBand,geometry,dynVars,ether,aux);
		
		for (i=0;i<ether.hilbertSize;i++) {
			//aux.eigM[i] += aux.eigOneBand[i];
		}
#ifdef MODEL_KONDO_DMS_MANYBANDS
		for (i=0;i<ether.numberOfOrbitals;i++) {
			temp = calcNumber(dynVars,geometry,ether,aux,i); // number of electrons for band i  
			s ="Number_Of_Band"+ttos(i)+"="+ttos(temp);
			io.historyPrint(s);
		}
#endif
		n_electrons=calcNumber (aux.eigOneBand,ether,aux);

	}
		
		
	s ="Number_Of_Electrons="+ttos(n_electrons);
	io.historyPrint(s);
		
	if (ether.tpem) temp=measure_energy (moment, ether,aux,tpemOptions);
	else temp=calcElectronicEnergy(aux.eigOneBand,ether,aux);
	
	s="Electronic Energy="+ttos(temp);
	io.historyPrint(s);
	temp2=calcSuperExchange(dynVars,geometry,ether);
	s="Superexchange="+ttos(temp2);
	io.historyPrint(s);
	temp += temp2;
	if (ether.isSet("tprime")) {
		temp2=directExchange2(dynVars,geometry,ether);
		s="Superexchange2="+ttos(temp2);
		io.historyPrint(s);
		temp += temp2;
	}
#ifdef MODEL_KONDO_INF_TWOBANDS		
	temp2 = calcPhononEnergy(dynVars,geometry,ether);
	s="PhononEnergy="+ttos(temp2);
	io.historyPrint(s);
	temp += temp2;
#endif		
	// total energy = electronic energy + superexchange + phonon energy
	s="TotalEnergy="+ttos(temp);
	io.historyPrint(s);
		
	if (ether.tpem) temp=measure_action (moment, ether,aux,tpemOptions);
	else temp=calcAction(aux.eigOneBand,ether,aux);
	s="Action="+ttos(temp);
	io.historyPrint(s);
	
	s="Number_Of_Holes="+ttos(ether.hilbertSize-n_electrons);
	io.historyPrint(s);

	s="ChemPot="+ttos(aux.varMu);
	io.historyPrint(s);
	
	if (ether.isSet("adjusttpembounds") || iter==0) {
		s="tpem_a="+ttos(aux.varTpem_a);
		io.historyPrint(s);
		s="tpem_b="+ttos(aux.varTpem_b);
		io.historyPrint(s);
	}
	
	for (i=0;i<n;i++) {
		tmp[i]=0.0;
		tmpweight[i]=0.0;
	}
	calcClasCor(geometry,dynVars,ether,tmp,tmpweight);
	if (iter==0 && ether.mpiRank==0 && ether.conc<ether.linSize) {
		cout<<"# weight = ";
		for (i=0;i<n;i++) {
			cout<<tmpweight[i]<<" ";
		}
		cout<<endl;
	}
	
	s="Cor[1]="+ttos(tmp[1]);
	io.historyPrint(s);

	if (geometry.isCubicType()) {
		int l = geometry.length();
		switch (d) {
			case 1:
				i=int(l*0.5);
			case 2:
				i+=int(l*l*0.5);
			case 3:
				i+=int(l*l*l*0.5);
		}
		temp=calcSq(tmp,geometry,i);
       
		s="Sq[pi]="+ttos(temp);
		io.historyPrint(s);
	
		if (d==2) {
			i=l*0.5;
			temp=calcSq(tmp,geometry,i);
			s="Sq[(pi,0)]="+ttos(temp);
			io.historyPrint(s);
        	}
	}
	for (i=0;i<n;i++) aux.clasCor[i] += tmp[i];
	
	if (iter==0) {
		s="Weight[0]="+ttos(tmpweight[0]);
		io.historyPrint(s);
	}
	temp= calcClasCor(dynVars,ether,0,4);
	s="CustomCorrelation="+ttos(temp);
        io.historyPrint(s);
	
	temp = calcMag(dynVars,ether);
	s="Mag2="+ttos(temp);
	io.historyPrint(s);
	
	kTpemMoments(moment,aux,ether);
	
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
	temp=calcPhononV(geometry,dynVars,ether);
	temp2=calcPhononV2(geometry,dynVars,ether);
	s="calcPhononV="+ttos(temp);
	io.historyPrint(s);
	s="calcPhononV2="+ttos(temp2);
	io.historyPrint(s);
	temp=temp2/ether.linSize - square(temp/ether.linSize);
	s="phononDeltaV="+ttos(temp);
	io.historyPrint(s);	
#endif		
	for (dof=0;dof<ether.classFieldList.size();dof++) {
		aux.nac[dof] = (double)aux.nac[dof]/ether.iterUnmeasured;
		s ="Accepted["+ttos(dof)+"]="+ttos(aux.nac[dof]);
		io.historyPrint(s);
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
	
	} else {
		temp=measure_kinetic(geometry,ether,aux,tpemOptions);
		if (ether.isSet("akw")) kTpemMomentsCl(geometry,ether,aux,tpemOptions);
		if (ether.isSet("optical")) kTpemMomentsOptical(geometry,ether,aux,tpemOptions);
		
	}


	s ="KineticEnergy="+ttos(temp);
	io.historyPrint(s);

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
		s="Conductance="+ttos(temp);
		io.historyPrint(s);
	}
	
	if (ether.bcsDelta0>0) {
		temp = calcBcsDelta2(dynVars,ether);
		s ="Delta2="+ttos(temp);
		io.historyPrint(s);
		calcBcsPhiCor(geometry,dynVars,ether,aux.bcsCorxx,0,0);
		calcBcsPhiCor(geometry,dynVars,ether,aux.bcsCorxy,0,1);
	}
	
	if (ether.isSet("saveall")) io.printSnapshot(dynVars,ether);
	if (ether.isSet("filelastconfig")) io.printSnapshot(dynVars,ether,1);

}				
