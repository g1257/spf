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




#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>

typedef  std::complex<double> myComplex;

struct My_params {
	std::string options,filename;
	double omegaMax,omegaStep,omegaMin,omega,beta,tpem_a,tpem_b,mu;
	int cutoff;
	std::vector<myComplex> coeffs;
};



// function that takes a double and a void * and returns a double for GSL library
typedef double  (*my_func_ptr)        (double,void *);

using namespace std;

void printValues(ostream &s,My_params const &params)
{
	s<<"options="<<params.options<<endl;
	s<<"filename="<<params.filename<<endl;
	s<<"beta="<<params.beta<<endl;
	s<<"tpem_a="<<params.tpem_a<<endl;
	s<<"tpem_b="<<params.tpem_b<<endl;
	s<<"mu="<<params.mu<<endl;
	s<<"cutoff="<<params.cutoff<<endl;
}
	
double core_integral(double a,double b,my_func_ptr my_function,My_params my_params)
{
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error; 

  
  gsl_function F;
  F.function = my_function;
  F.params = &my_params;

  gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,
                        w, &result, &error); 

  gsl_integration_workspace_free(w);
  return result;
}

double fermi(double x)
{
	double res;
	if (x>50) return 0;
	if (x<-50) return 1;
	res = 1.0/(1.0+exp(x));
	return res;
}

double tpem_chebyshev(int m,double x) {
	double tmp;
	int p;
	if (m==0) return 1;
	if (m==1) return x;
	
	if ((m%2)==0) {
		p=m/2;
		tmp=tpem_chebyshev(p,x);
		return (2*tmp*tmp-1);
	}
	else {
		p=(m-1)/2;
		return (2*tpem_chebyshev(p,x)*tpem_chebyshev(p+1,x)-x);
	}
}

double sigma_integrand(double x,void *params)
{
	My_params *p = (My_params *)params;
	myComplex tmp = myComplex(0.0,0.0);
	double omega = p->omega;
	int cutoff = p->cutoff;
	double beta = p->beta;
	int l,m;
	double omega_norm,x_norm;
	
	for (l=0;l<cutoff;l++) {
		for (m=0;m<cutoff;m++) {
			x_norm = (x-p->tpem_b)/p->tpem_a;
			omega_norm = (omega - p->tpem_b)/p->tpem_a;
			tmp += p->coeffs[l+m*cutoff] * tpem_chebyshev(m,x_norm) * tpem_chebyshev(l,x_norm+omega_norm);
		}
	} 
	tmp *= (fermi(beta * (x-p->mu) )-fermi(beta * (x+omega-p->mu) ))/omega;
	return real(tmp);
}

	
	
void loadCoeffs(My_params &params)
{
	int i,index;
	double tmp,tmp2;
	
	
	ifstream fin(params.filename.c_str());
	if (!fin || fin.bad()) {
		cerr<<"Cannot open or read file "<<params.filename<<endl;
		exit(1);
	}
	fin>>params.cutoff;
	fin>>params.beta;
	fin>>params.mu;
	fin>>params.tpem_a;
	fin>>params.tpem_b;
	
	int cutoff = params.cutoff;
	
	cerr<<(cutoff*cutoff)<<" <== CUTOFF\n";
	
	for (i=0;i<cutoff*cutoff;i++) {
		fin>>index;
		fin>>tmp;
		fin>>tmp2;
		if (i!=index) {
			cerr<<"Expected "<<i<<" but got "<<index<<endl;
			exit(1);
		}
		params.coeffs.push_back(myComplex(tmp,tmp2));
	}
	fin.close();
}	
	
void initValues(char *argv[],My_params &params)
{	
	
	ifstream fin(argv[1]);
	if (!fin || fin.bad()) {
		cerr<<"Cannot open or read file "<<argv[1]<<endl;
		exit(1);
	}
	fin>>params.options;
	fin>>params.filename;
	fin>>params.omegaMin;
	fin>>params.omegaMax;
	fin>>params.omegaStep;
	
	
	fin.close();
	
}

void my_handler (const char * reason, const char * file, int line, int gsl_errno
)
{
        fprintf(stderr,"WARNING (Ignore), %s, file=%s, line=%d, code=%d\n",
                        reason,file,line,gsl_errno);
}

int main(int argc,char *argv[])
{
	My_params params;
	double omega,tmp;
	
	gsl_set_error_handler (&my_handler);
		
	initValues(argv,params);
	loadCoeffs(params);
	printValues(cerr,params);
	
	for (omega=params.omegaMin;omega<params.omegaMax;omega+=params.omegaStep) {
		params.omega = omega;
		if (omega>=2*params.tpem_a) break;
		tmp = core_integral(params.tpem_b- params.tpem_a,params.tpem_b + params.tpem_a -omega,sigma_integrand,params);
		cout<<omega<<" "<<tmp<<endl;
	}
	return 0;
}
