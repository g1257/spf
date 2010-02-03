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




#include "common.h"
#include "powell.h"

extern double calcSuperExchange(DynVars const &dynVars, Geometry const &geometry, Parameters const &ether);
 
Geometry sgeometry;
DynVars sdynVars;
Aux saux;
Parameters sether;

using std::cout;
	
double gsEnergyInternal(vector<double> &Theta,Geometry const &geometry,DynVars &dynVars,Aux &aux,Parameters const &ether)
{
	int i;
	double s=0.;
	
	for (i=0;i<geometry.volume();i++) {
		/* if (cos(Theta[i])>0) dynVars.theta[i]=0.0;
		else  dynVars.theta[i]=M_PI; */
		dynVars.theta[i]=Theta[i];
		dynVars.phi[i]=0.0;
	} 
	
	
		
	setupHamiltonian(aux.matrix,geometry,dynVars,ether,aux,0);
	diag(aux.matrix,aux.eigOneBand,'n');
	for (i=0;i<ether.carriers;i++) {
		s += aux.eigOneBand[i];
	}
	s+= calcSuperExchange(dynVars,geometry, ether);
	return s;
}

double gsEnergy(vector<double> &theta_nonzero)
{
	vector<double> theta;
	unsigned int j=0;
	for (int i=0;i<sgeometry.volume();i++) {
		theta.push_back(0.0);
		if (j>=theta_nonzero.size()) {
			cerr<<"Problem at "<<__FILE__<<" "<<__LINE__<<endl;
			exit(1);
		}
		if (sether.modulus[i]!=0) theta[i]=theta_nonzero[j++];
	}
	return gsEnergyInternal(theta,sgeometry,sdynVars,saux,sether);
}
	
void calcGroundState(Geometry const &geometry,DynVars &dynVars,Parameters const &ether,Aux &aux)
{
	size_t n=ether.conc;
	vector<double> p(n);
	vector<vector<double> > xi;
	vector<double> vecTmp(n);
	double ftol = 1e-4;
	int iter=0;
	double fret=0;
	vector<double> theta_nonzero;
	
	for (int i=0;i<geometry.volume();i++) {
		if (ether.modulus[i]!=0) theta_nonzero.push_back(dynVars.theta[i]);
	}
	
	if (theta_nonzero.size()!=n) {
		cerr<<"Problem at "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			vecTmp[j]=0.0; 
			if (i==j) vecTmp[j]=1;
		}		
		xi.push_back(vecTmp); // versors
		p[i]=theta_nonzero[i];
	}
	
	
	sgeometry= geometry;
	sdynVars=dynVars;
	saux=aux;
	sether=ether;
	cerr<<"In powell with n="<<n<<"\n";
	powell(p,xi,n,ftol,iter,fret, gsEnergy);
	vectorPrint(p,"gs",cerr);
	for (int i=0,j=0;i<geometry.volume();i++) {
		dynVars.theta[i]=0.0;
		dynVars.phi[i]=0.0;
		if (size_t(j)>=p.size()) {
			cerr<<"Problem at "<<__FILE__<<" "<<__LINE__<<endl;
			exit(1);
		}
		if (ether.modulus[i]!=0) dynVars.theta[i]=p[j++];
	}
	
	cout<<"iter="<<iter<<" and fret="<<fret<<endl; 
}


