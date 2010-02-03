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
#include "lanczosplus.h"
#include "ridders.h"

/********** ADJUST BOUNDS OF TPEM SPECTRUM USING LANCZOS ***********************/
void tpem_sparse_scale(tpem_sparse *matrix,double a,double b)
{
	
	for(size_t i=0;i<matrix->rank; i++) {
	  for(size_t j= (matrix->rowptr[i]); j < (matrix->rowptr[i+1]); j++) {
	    if ((matrix->colind[j]) == i)
	      matrix->values[j] = (matrix->values[j]-b)/a;
	    else 
	     matrix->values[j] = matrix->values[j]/a;
		
	  }
	}
}
	
int tpemAdjustBounds(tpem_sparse *matrix,Parameters const &ether, Aux &aux)
{	
	SparseMatrix mat (matrix);
	vector<MatType> z;
	
	int max_step=100; //hardcoded
	
	
 	lanczos (0, mat, max_step, aux.varEnergy2, z);
	lanczos (1, mat, max_step, aux.varEnergy1, z);
	
	double temp = aux.varEnergy2 - aux.varEnergy1;
	
	if (temp < 0) 
	  { 
	    cerr << "E_max < E_min !!!!" << endl;
	    return 1;
	  }
	aux.varEnergy1 -= 0.05*temp; 
	aux.varEnergy2 += 0.05*temp;
	
	
	
	aux.varTpem_a = 0.5 * (aux.varEnergy2 - aux.varEnergy1);
	aux.varTpem_b = 0.5 * (aux.varEnergy2 + aux.varEnergy1);
	
	tpem_sparse_scale(matrix,aux.varTpem_a,aux.varTpem_b);
	
	return 0;
}

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/



double nOfElectrons(double mu,double beta,vector<double> const &eig)
{
	unsigned int i;
	double n_electrons;
	
	n_electrons=0;
	for (i=0;i<eig.size();i++) {
		n_electrons += fermi((eig[i]-mu)*beta);
		//cerr<<"Here is i="<<i<<" and eig="<<eig[i]<<"and factor="<<(mu-eig[i])*beta<<"\n";
	}
	return n_electrons;
}	

// Derivative of n(mu) with respect to mu
double nOfElectPrime(double mu,double beta,vector<double> const &eig)
{
	unsigned int i;
	double prime;
	
	prime=0;
	for (i=0;i<eig.size();i++) {
		prime -= fermiPrime((eig[i]-mu)*beta)*beta;
	}
	return prime;
}


void adjChemPot(vector<double> const &eig,Parameters const &ether,Aux &aux)
{
	double denom,mu=aux.varMu,tmp;
	bool converged=false;
	int iter;
	int maxIter=100,n0=ether.carriers;
	double tolerance=1.0e-3;
	
	for (iter=0;iter<maxIter;iter++) {
		denom=nOfElectPrime(mu,ether.beta,eig);
		tmp = nOfElectrons(mu,ether.beta,eig);
		if (fabs(denom)<1e-3) break;
		mu = mu -(tmp-n0)/denom;
		//cerr<<"electrons="<<tmp<<" n0="<<n0<<endl;
         //       cerr<<"Denom="<<denom<<" Proposed mu="<<mu<<endl;
                
		if (fabs(nOfElectrons(mu,ether.beta,eig)-n0)<tolerance) {
			converged=true;
			break;
		}
	}
	//if (!converged) cerr<<"Failed to converged after "<<iter<<" iterations.\n";
	if (converged) aux.varMu=mu;
}

// ADJUST CHEMICAL POTENTIAL FOR TPEM ****************/

double nElectronsTpem (double x)
{
	double	ret;
	double	a, b, mu, beta;
	
	tmpValues(a,b,mu,beta,1); 
    	ret = 0.5 * (1.0 - tanh (0.5 * beta * (a * x + b - mu)));
	return ret;
}

/*
double nElectronsPrimeTpem (double x)
{
	double	ret;
	double	a, b, mu, beta;
		
	tmpValues(a,b,mu,beta,1);
	//	cerr << a << b << mu << beta << endl;
	ret = -beta*0.25*(1-square(tanh (0.5 * beta * (a * x + b - mu))));
	return ret;
}*/

double measure_nOfElectrons (double mu,vector<double> const &moment,double beta,double a,double b, TpemOptions const
	&tpemOptions)
{
	int cutoff=moment.size();
	vector<double> coeffs(cutoff);
 	
	tmpValues(a,b,mu,beta,0);
	tpem_calculate_coeffs (coeffs, nElectronsTpem,tpemOptions);
	return tpem_expansion (moment, coeffs);

}

/*
double measure_nOfElectPrime (double mu,vector<double> const &moment,
	Parameters const &ether,Aux &aux)
{
	int cutoff = ether.tpem_cutoff;
	vector<double> coeffs(cutoff);
	double a=aux.varTpem_a,b=aux.varTpem_b,beta=ether.beta;
   	TpemOptions tpemOptions;
	
	tpemOptionsFill(tpemOptions,ether);
	tmpValues(a,b,mu,beta,0);
	tpem_calculate_coeffs (coeffs, nElectronsPrimeTpem,tpemOptions);
	return  tpem_expansion (moment, coeffs);

  
}
*/

struct My_params
{
	vector<double> moments;
	double a,b,beta;
	int cutoff,carriers;
	TpemOptions tpemOptions;
};

// n(mu)-n0
double fForAdjMu(double x,My_params const &params)
{
	int i,cutoff = params.cutoff;
	double ret;
	vector<double> moment(cutoff);
	
	for (i=0;i<cutoff;i++) moment[i]=params.moments[i];
	
	ret = measure_nOfElectrons (x,moment,params.beta,params.a,params.b,params.tpemOptions);
	//cerr<<"Evaluating at x="<<x<<" is ret="<<ret<<endl;
	ret -= params.carriers;
	return ret;
	
}

void adjChemPotTpem(Parameters const &ether, Aux &aux,TpemOptions const
&tpemOptions)
{
	double mu,lowband=ether.energy1+0.1;
	int info=1;
	My_params params;
	int cutoff=ether.tpem_cutoff;
	vector<double> moment(cutoff);
	int iter;
      
	if (ether.tpem==1)  {
		for (iter=0;iter<cutoff;iter++) moment[iter] = aux.curMoments[iter]; 
	} else {			
		tpem_calculate_moment(aux.sparseMatrix[0],moment,tpemOptions);
	}
      //vectorPrint(moment,"momentsaccurate",std::cout);
      //vectorPrint(aux.curMoments,"momentsfast",std::cout);
	
	for (iter=0;iter<cutoff;iter++) params.moments.push_back(moment[iter]);
	
	
	params.a=aux.varTpem_a;
	params.b=aux.varTpem_b;
	params.beta=ether.beta;
	params.cutoff=ether.tpem_cutoff;
	params.carriers=ether.carriers;
	params.tpemOptions=tpemOptions;
	
 	if (tpemOptions.tpem_rank==0) mu=root_main(fForAdjMu,lowband,0,params,info); // from ridders.h
	tpem_bcast(&mu,tpemOptions);
	tpem_bcast(&info,tpemOptions);
	if (info==0) {
		//std::cerr<<"ether.mpiRank="<<ether.mpiRank<<" "<<info<<" "<<mu<<endl;		
		aux.varMu=mu;
		
	}
	
}


