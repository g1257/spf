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




#include "ridders.h"

using namespace std;


double maxElement(vector<double> const &v);
double minElement(vector<double> const &v);

int root_count(vector<double> const &ff)
{
	int c=0,i,n=ff.size();
	
	for (i=0;i<n;i++) if (ff[i]>0) c++;
	return c;
}

double root_aux_find(vector<double> const &ff)
{

	int c=root_count(ff);
		
	// c==1 then discard negative, if cc==2 then discard positive
	if (c==1) return minElement(ff);
	if (c==2) return maxElement(ff);
	cerr<<"c="<<c<<" in root_aux_find\n";
	exit(1);
	return 0; // to avoid compiler warning
}	
	


void root_aux_discard(vector<double> &a,vector<double> const &ff,double v)
{
	vector<double> b(3);
	int i,k=0,n=a.size();
	
	for (i=0;i<n;i++) b[i]=0;
	
	for (i=0;i<n;i++) {
		if (ff[i]!=v) {
			b[k]=a[i];
			k++;
		}
	}
	for (i=0;i<n;i++) a[i]=b[i];
}

double root_main(my_fftc_ptr f,double a1,double a2,My_params const &my_params,int &info)
{
	vector<double> a(3);
	vector<double> ff(3);
	int i,max_counter=200;
	double v,tolerance=1e-3;
	int counter;
	a[0]=a1; a[1]=a2; 
	info=1;
	
	counter=0;
	//cerr<<"Testing f(0)"<<f(0,my_params)<<" and f(-8)="<<f(-8.,my_params)<<endl;
	while(counter<max_counter) {
		a[2]=(a[0]+a[1])*0.5; // medium point
		//cerr<<"Studying interval "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
		for (i=0;i<3;i++) {
			ff[i]=f(a[i],my_params); // y-axis points
			if (fabs(ff[i])<tolerance) {
				info=0;
				return a[i];
			}
		}
		//cerr<<"Studying interval with f= "<<ff[0]<<" "<<ff[1]<<" "<<ff[2]<<endl;
		if (root_count(ff)<1 || root_count(ff)>2) return 0;
		v=root_aux_find(ff);
		//cerr<<"root to discard is v="<<v<<endl;
		root_aux_discard(a,ff,v); // discard it
		counter++;
	}	
	if (counter>=max_counter) {	 
		cerr<<"Convergence NOT achieved EVEN after "<<counter<<" steps.\n";
	}

	return 0;
}

