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




#include "powell.h"

#define TINY 1.0e-25
#define ITMAX 200 
#define TOL 2.0e-4
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GLIMIT 100.0 
#define GOLD 1.618034 

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

int ncom;
double (*nrfunc)(std::vector<double> &);
std::vector<double> pcom,xicom;

using namespace std;

void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, double (*func)(double))  
{ 
	double ulim,u,r,q,fu,dum; 
	fa=(*func)(ax); 
	fb=(*func)(bx); 
	if (fb > fa) { 
		SHFT(dum,ax,bx,dum)
		SHFT(dum,fb,fa,dum) 
	} 
	cx=bx+GOLD*(bx-ax); 
	fc=(*func)(cx); 
	while (fb > fc) { 
		r=(bx-ax)*(fb-fc); 
		q=(bx-cx)*(fb-fa);
		u=(bx)-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(bx)+GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) { 
			fu=(*func)(u); if (fu < fc) { 
			ax=(bx); 
			bx=u; 
			fa=(fb); 
			fb=fu; 
			return; 
		} else if (fu > fb) { 
			cx=u; fc=fu; 
			return; 
		} 
		u=(cx)+GOLD*(cx-bx); 
		fu=(*func)(u); 
		} else if ((cx-u)*(u-ulim) > 0.0) { 
			fu=(*func)(u); 
			if (fu < fc) { 
				SHFT(bx,cx,u,cx+GOLD*(cx-bx)) 
				SHFT(fb,fc,fu,(*func)(u)) 
			} 
		} else if ((u-ulim)*(ulim-cx) >= 0.0) { 
			u=ulim; 
			fu=(*func)(u); 
		} else { 
			u=(cx)+GOLD*(cx-bx); 
			fu=(*func)(u); 
		}	 
		SHFT(ax,bx,cx,u) 
		SHFT(fa,fb,fc,fu) 
	} 
} 

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double &xmin) 
{
int iter;
double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
double e=0.0;
 
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx; 
	fw=fv=fx=(*f)(x);
 
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS); 
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) { 
			xmin=x; 
			return fx;
		}
		if (fabs(e) > tol1) { 
			r=(x-w)*(fx-fv); 
			q=(x-v)*(fx-fw); 
			p=(x-v)*q-(x-w)*r; 
			q=2.0*(q-r); 
			if (q > 0.0) p = -p; 
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
			else { 
				d=p/q; 
				u=x+d; 
				if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x); 
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
		} 
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d)); fu=(*f)(u); 
		if (fu <= fx) { 
			if (u >= x) a=x; 
			else b=x; 
			SHFT(v,w,x,u) 
			SHFT(fv,fw,fx,fu) 
		} else { 
			if (u < x) a=u; 
			else b=u; 
			if (fu <= fw || w == x) { 
				v=w;
				w=u;
				fv=fw;
				fw=fu; 
			} else if (fu <= fv || v == x || v == w) { 
				v=u;
				fv=fu; 
			} 
		} 
	} 
	cerr<<"Too many iterations in brent"<<endl;
	exit(1); 
	xmin=x; 
	return fx; 
} 
 
double f1dim(double x)
{
	int j;
	double f;
	vector<double> xt(ncom);
	
	
	for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	
	return f;
}


void linmin(vector<double> &p , vector<double> &xi, int n, double &fret, double (*func)(vector<double> & ))
{
	
	int j;
	double xx,xmin,fx,fb,fa,bx,ax; 
	
	ncom=n;
	pcom.insert(pcom.begin(),n,0.0);
	xicom.insert(xicom.begin(),n,0.0); 
	
	nrfunc=func;
	for (j=0;j<n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(ax,xx,bx,fa,fx,fb,f1dim); //FIXME
	fret=brent(ax,xx,bx,f1dim,TOL,xmin); // FIXME
	for (j=0;j<n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	pcom.clear();
	xicom.clear(); 
} 
 
void powell(vector<double> &p, vector<vector<double> > &xi, int n, double ftol, int &iter, 
double  &fret, double (*func)(vector<double> &)) 
/* Minimizationof a function func of n variables. 
Input consists of an initial starting point p[1..n]; an initial matrix xi[1..n][1..n], 
whose columns contain the initial set of directions (usually the nunit vectors); and ftol, 
the fractional tolerance in the function value such that failure to decrease by more than this amount 
on one iteration signals doneness. On output, p is set to the best point found, 
xi is the then-current direction set, 
fret is the returned function value at p,and iter is the number of iterations taken. The routine linmin 
is used. 
*/
{ 
	int i,ibig,j; 
	double del,fp,fptt,t; 
	vector<double> pt(n),ptt(n),xit(n);
	
	/* pcom.insert(pcom.begin(),n,0.0);
	xicom.insert(xicom.begin(),n,0.0); */
	
	fret=(*func)(p);
	for (j=0;j<n;j++) pt[j]=p[j]; 
	/* Savetheinitial point. */
	for (iter=1;;++iter) {
		cerr<<"powell: iteration "<<iter<<endl;
		fp=(fret); ibig=0;
		del=0.0;
		/* Will bethebiggest functiondecrease. */
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i]; 
				fptt=(fret);
				linmin(p,xit,n,fret,func);
				if (fptt-(fret) > del) {
					del=fptt-(fret);
					ibig=i;
			}
		}
		if (2.0*(fp-(fret)) <= ftol*(fabs(fp)+fabs((double)fret))+TINY) {
			xit.clear();
			ptt.clear();
			pt.clear();
			cerr<<"returning...fp="<<fp<<"\n";
			return;
		}
		if (iter == ITMAX) {
			cerr<<"powell exceeding maximum iterations.";
			exit(1);
		}
		for (j=0;j<n;j++) { 
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(fret)+fptt)*sqrt(fp-(fret)-del)-del*sqrt(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=0;j<n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	} 
} 
