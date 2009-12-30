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




#include "parameters.h"

using namespace std;

void Parameters::print(ostream &s,char prefix)
{
	int i,j;
	
#ifdef USE_MPI
	if (mpiTpemRank==0) typeofmodel = typeofmodel + string("+MPI");
#endif	
	s<<prefix<<"This file was created with the SPF program\n";
	s<<prefix<<"VERSION="<<version<<endl;
	s<<prefix<<"TYPE_OF_MODEL="<<typeofmodel<<endl;
	s<<prefix<<"COMPILATION_DATE="<<__DATE__<<endl;
	s<<prefix<<"options="<<options<<endl;
	s<<prefix<<"lattice="<<typeoflattice<<endl;
	s<<prefix<<"iterTherm="<<iterTherm<<endl;
	s<<prefix<<"iterEffective="<<iterEffective<<endl;
	s<<prefix<<"iterUnmeasured="<<iterUnmeasured<<endl;
	s<<prefix<<"linSize="<<linSize<<endl;
	if (carriers>0) s<<prefix<<"carriers="<<carriers<<endl;
	s<<prefix<<"beta="<<beta<<endl;
	s<<prefix<<"mcflag="<<mcflag<<endl;
	s<<prefix<<"window="<<window<<endl;
	s<<prefix<<"Rootname="<<rootname<<endl;
	s<<prefix<<"startType="<<startType<<endl;
	if (startType==4 || startType==6) {
		s<<prefix<<"dynVarInputFile="<<dynVarsInputFile<<endl;
	}
	if (startType==6) {
		s<<prefix<<"startLevel="<<startLevel<<endl;
	}
	s<<prefix<<"HistogramSteps="<<histSteps<<endl;
	s<<prefix<<"HistgramBounds="<<hist1<<" "<<hist2<<endl;
	s<<prefix<<"conc="<<conc<<endl;
	s<<prefix<<"hoppings=";
	for (i=0;i<D;i++) s<<" "<<hoppings[i];
	s<<endl;
#ifdef MODEL_KONDO_FINITE	
	s<<prefix<<"Couplings="<<J[0]<<endl;
#endif
#ifdef MODEL_KONDO_DMS_THREEBANDS
	s<<prefix<<"Hoppings=";
	for (i=0;i<9*D;i++) {
		s<<bandHoppings[i]<<" ";
	}
	s<<endl; 
	/* s<<prefix<<"Couplings=";
	for (i=0;i<2;i++) {
		s<<J[i]<<" ";
	}
	s<<endl; */
#endif
#ifdef MODEL_KONDO_DMS_MANYBANDS
	s<<prefix<<"Hoppings=";
	for (i=0;i<bandHoppings.size();i++) {
		s<<bandHoppings[i]<<" ";
	}
	s<<endl;
	s<<prefix<<"Couplings=";
	for (i=0;i<J.size();i++) {
		s<<J[i]<<" ";
	}
	s<<endl;
#endif
#ifdef MODEL_KONDO_INF_ONEBAND
	s<<prefix<<"Couplings=infinity"<<endl;	
#endif
#ifdef MODEL_KONDO_INF_TWOBANDS
	s<<prefix<<"Couplings=infinity"<<endl;	
#endif	
	s<<prefix<<"JAF coupling(s)=";
	j=JafVector.size();
	if (j>8 && !isSet("verbose")) j=8;
	if (isSet("jafvector")) {
		s<<jafFile;
		for (i=0;i<j;i++) s<<JafVector[i]<<" ";
	} else {
		for (i=0;i<j;i++) s<<JafVector[i]<<" "; 
	}
	
	s<<endl;
	if (isSet("havepotential")) {
		
		s<<prefix<<"PotentialFile="<<potentialFile<<" ";
	} else {
		s<<prefix<<"Potential ";
	}
	
	j=potential.size();
	if (j>8 && !isSet("verbose")) j=8;
	for (i=0;i<j;i++) s<<potential[i]<<" ";
	s<<endl;
		
	
	
	if (isSet("magneticfield")) {
		s<<prefix<<"Magnetic Field=";
		for (i=0;i<3;i++) s<<magnetic[i]<<" ";
		s<<endl;
	}
	
	if (isSet("tprime")) {
		s<<prefix<<"TPRIME="<<tprime<<endl;
		s<<prefix<<"JAFPRIME="<<jafprime<<endl;
	}
	if (isSet("tsecond")) s<<prefix<<"TSECOND="<<tsecond<<endl;

	s<<prefix<<"ENERGY_LOWER_BOUND="<<energy1<<endl;
	s<<prefix<<"ENERGY_UPPER_BOUND="<<energy2<<endl;
	if (tpem) {
		s<<prefix<<"TPEM_CUTOFF="<<tpem_cutoff<<endl;
		s<<prefix<<"TPEM_EPSPROD="<<tpem_epsProd<<endl;
		s<<prefix<<"TPEM_EPSTRACE="<<tpem_epsTrace<<endl;
		
	}		
#ifdef MODEL_KONDO_INF_TWOBANDS
	s<<prefix<<"#BandHoppings ";
	for (i=0;i<bandHoppings.size();i++) {
		s<<bandHoppings[i]<<" ";
	}
	s<<endl;
	s<<prefix<<"#EJTCouplings ";
	for (i=0;i<phononEjt.size();i++) {
		s<<phononEjt[i]<<" ";
	}
	s<<endl;
	s<<prefix<<"#EdCouplings ";
	for (i=0;i<phononEd.size();i++) {
		s<<phononEd[i]<<" ";
	}
	s<<endl;
	s<<prefix<<"MaxPhonons "<<maxPhonons<<endl;
#endif	
#ifdef MODEL_KONDO_PHONONS
	s<<prefix<<"Lambda="<<phononLambda<<endl;
#endif
#ifdef MODEL_KONDO_PHONONS_EX	
	s<<prefix<<"Alpha="<<phononAlpha<<"\n";
	s<<prefix<<"Delta="<<phononDelta<<"\n";
	s<<prefix<<"GJ="<<phononGj<<"\n";
#endif
	if (bcsDelta0>0) {
		s<<prefix<<"BCSDELTA0="<<bcsDelta0<<endl;
		s<<prefix<<"BCSV=FIXME (sorry)"<<endl;
	}
	s<<prefix<<"Modulus=";
	if (conc<linSize) {
		for (i=0;i<linSize;i++) if (modulus[i]>0) s<<" "<<i;
	} else {
		s<<"Concentrated";
	}
	s<<endl;
#ifdef USE_MPI	
	s<<prefix<<"NumberOfProcessors="<<mpiNop1<<" times "<<mpiNop2<<"\n";
#endif	
	
}


Parameters::Parameters()
{
	iterTherm=iterEffective=iterUnmeasured=0;
	window=beta=0;
	linSize=mcflag=startType=conc=0;
	mpiRank=0;
	mpiNop1=mpiNop2=1;
	version="v6.4. Released: March 29, 2006.";
}

bool Parameters::isSet(char const *what) const
{
	std::string s1(options);
	std::string s2(what);
	if (s1.find(s2)==std::string::npos) return false;
	return true;
}

void Parameters::setOption(char const *what) 
{
	if (isSet(what)) {
		cerr<<"Option "<<what<<" already set\n";
		return;
	}
	options = options + what;
}


void Parameters::welcome()
{
	if (mpiRank>0) return;
	string ansiThis = ansiYellow;
	string versionColor = "yellow";
	std::cout<<"\n\n************************************\n";
	if (!isSet("noansicolors")) std::cout<<ansiThis;
	std::cout<<"SPF "<<version<<" ("<<versionColor<<")\n";
	if (!isSet("noansicolors")) std::cout<<ansiReset;
	std::cout<<"TYPE_OF_MODEL="<<typeofmodel<<endl;
	std::cout<<"VERSION="<<version<<" Compiled:"<<__DATE__<<endl;
	std::cout<<"For more info see: "
	"http://mri-fre.ornl.gov/spf\n";
	std::cout<<"\n";
}

