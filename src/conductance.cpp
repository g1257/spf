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




#include "conductance.h"

using namespace std;

#ifdef MODEL_KONDO_FINITE
#define MODEL_KONDO_2N
#endif
#ifdef MODEL_KONDO_INF_TWOBANDS
#define MODEL_KONDO_2N
#endif
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
#define MODEL_KONDO_1N
#endif
#ifdef MODEL_KONDO_INF_ONEBAND
#define MODEL_KONDO_1N
#endif


complex<double> **matrix_alloc (unsigned int row, unsigned int col)
{
        size_t i;
        complex<double> **this1a = new complex<double>*[row];

        *this1a = new complex<double>[row * col];
        for (i = 1; i < row; i++)
                this1a[i] = this1a[i - 1] + col;
        return this1a;
}

void matrix_free (complex<double> **this1a)
{
        if (this1a && *this1a) delete [] *this1a, *this1a = NULL;
        if (this1a) delete []  this1a, this1a = NULL;
}

#ifndef MODEL_KONDO_2N
#ifndef MODEL_KONDO_1N
double conductance3d(MyMatrix<MatType> const & matrix,double mu, int n,int dim,int l1,
int maxiter,double maxerror)
{
	static int firstcall=1;
	if (firstcall) {
		std::cerr<<"WARNING: conductance is not implemented for this model and will always return -1.\n";
		firstcall=0;
	}
	return -1; // not implemented
}
#endif
#endif

#ifdef MODEL_KONDO_2N
double conductance3d(MyMatrix<MatType> const & matrix,double mu, int n,int dim,int l1,
int maxiter,double maxerror)
{
	static	int	firstcall = 1;
	static int volume, *pivot;
	static	Complex	**h0, **h1, **t0, **t1, **x;
	static	Complex	**sigma_right, **sigma_left, **green_plus;
	char c123 = 'n';
	int	 m, i, j, x0, x1, y0, y1, z0, iter;
	double	eta = 0.0, a,b, error, cond;
	Complex	one = Complex (1.0, 0.0), zero = Complex (0.0, 0.0), e;
	
	if (n == 0) {
		if (pivot) delete [] pivot, pivot = NULL;
		matrix_free (h0), h0 = NULL;
		matrix_free (h1), h1 = NULL;
		matrix_free (t0), t0 = NULL;
		matrix_free (t1), t1 = NULL;
		matrix_free (x), x = NULL;
		matrix_free (sigma_right), sigma_right = NULL;
		matrix_free (sigma_left),  sigma_left  = NULL;
		matrix_free (green_plus),  green_plus  = NULL;
		firstcall = 1;
		return -1;
	}
	
	/* matrix = (tpem_sparse *) this->data; */
	m = n * n;
	if (firstcall) {
		firstcall = 0;
		/* dimension = *(int *) parser_get ("dimension");
		volume    = *(int *) parser_get ("volume"); */
		pivot = new int [n]; 
		volume=n/2;
		
		h0 = matrix_alloc (n, n);
		h1 = matrix_alloc (n, n);
		t0 = matrix_alloc(n, n);
		t1 = matrix_alloc(n, n);
		x  = matrix_alloc(n, n);
		sigma_right = matrix_alloc (n, n);
		sigma_left  = matrix_alloc(n, n);
		green_plus  = matrix_alloc(n, n);
	}
	
	/* intralayer (resp. interlayer) hamiltonian h0 (resp. h1) */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			x[i][j] = zero;
			h0[i][j]=h1[i][j]=0;
		}
	}
	
	
	
	switch (dim) {
	case 3:
		for (i = 0; i < l1; i++) {
		for (j = 0; j < l1; j++) {
			for (int k=0;k<l1;k++) {
				x0=k*l1*l1+j*l1+i;
				x[x0][x0]=x[x0+volume][x0+volume]=k;
			}
		}
	}
	
	for (i = 0; i < l1; i++) {
	for (j = 0; j < l1; j++) {
	for (int k = 0; k < l1; k++) {
		x1 = k*l1*l1 + j*l1 +i;
		for (x0=0;x0<l1;x0++) {
			for (y0=0;y0<l1;y0++) {
				for (z0=0;z0<l1;z0++) {
					y1=z0*l1*l1+y0*l1+x0;
					if (z0==0 && k==l1-1) {
						h1[x1][y1]=matrix(x1,y1);
						h1[x1+volume][y1+volume]=matrix(x1+volume,y1+volume);
						h1[x1+volume][y1]=matrix(x1+volume,y1);
						h1[x1][y1+volume]=matrix(x1,y1+volume);
					}
					else if (k!=0 || z0!=l1-1) {
						h0[x1][y1]=matrix(x1,y1);
						h0[x1+volume][y1+volume]=matrix(x1+volume,y1+volume);
						h0[x1+volume][y1]=matrix(x1+volume,y1);
						h0[x1][y1+volume]=matrix(x1,y1+volume);
						
					}
				}
			}
		}
	}}}
	break;
	case 2:
		for (i = 0; i < l1; i++) {
		for (j = 0; j < l1; j++) {
			x0=j*l1+i;
			x[x0][x0]=x[x0+volume][x0+volume]=j;
		}
		}
	
	for (i = 0; i < l1; i++) {
	for (j = 0; j < l1; j++) {
		x1 =j*l1 +i;
		for (x0=0;x0<l1;x0++) {
			for (y0=0;y0<l1;y0++) {
					y1=y0*l1+x0;
					if (y0==0 && j==l1-1) {
						h1[x1][y1]=matrix(x1,y1);
						h1[x1+volume][y1+volume]=matrix(x1+volume,y1+volume);
						h1[x1+volume][y1]=matrix(x1+volume,y1);
						h1[x1][y1+volume]=matrix(x1,y1+volume);
					}
					else if (j!=0 || y0!=l1-1) {
						h0[x1][y1]=matrix(x1,y1);
						h0[x1+volume][y1+volume]=matrix(x1+volume,y1+volume);
						h0[x1+volume][y1]=matrix(x1+volume,y1);
						h0[x1][y1+volume]=matrix(x1,y1+volume);
						
					}
				}
			}
	}}
	break;
	case 1:
	cerr<<"Cannot calculate conductance in 1d\n";
	return -1;
	break;
}
	
#ifndef _AIX
	/* velocity operator */
	
	zgemm_ (&c123, &c123, &n, &n, &n, &one, *h0, &n, *x, &n, &zero, *t0, &n);
    zgemm_ (&c123, &c123, &n, &n, &n, &one, *x, &n, *h0, &n, &zero, *t1, &n);
#else
	/* zgemm(*h0,&n,"n",*x,&n,"c",*t0,&n,&n,&n,&n,NULL,0);
	zgemm(*x,&n,"n",*h0,&n,"n",*t1,&n,&n,&n,&n,NULL,0); */
	zgemm ('n', 'n', &n, &n, &n, &one, *h0, &n, *x, &n, &zero, *t0, &n);
    zgemm ("n", "n", &n, &n, &n, &one, *x, &n, *h0, &n, &zero, *t1, &n);
#endif
	
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
        		x[i][j] = Complex((double)(imag(t0[i][j] - t1[i][j])) /
				(l1 - 1.0),
		(double)(real(t1[i][j] - t0[i][j])) /
				(l1 - 1.0));
		}
	
	e = Complex (mu, eta);

	/* right self-energy */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
          		sigma_right[i][j] = zero;
		sigma_right[i][i] = Complex (0.0, -1.0);
	}
	
	error = 1.0;
	for (iter = 0; iter < maxiter && error > maxerror; iter++) {

		for (i = 0; i < n; i++) 
			for (j = 0; j < n; j++) {
				if (i!=j)
          			t1[i][j] =  0.0- h0[i][j]- sigma_right[i][j];
				else
					t1[i][j] = e - h0[i][j]- sigma_right[i][j];
          			
			}
#ifndef _AIX	
		zgetrf_ (&n, &n, *t1, &n, pivot, &x0);
		zgetri_ (&n, *t1, &n, pivot, *t0, &m, &x0);
		/* t1 = h1 * t1 * h1^H */
		zgemm_ (&c123,&c123, &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm_ (&c123,&c123, &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
	
#else
		for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) green_plus[i][j]=zero;
			green_plus[i][i]=one; 
		}
		zgef(*t1,&n,&n,pivot);
		zgesm("n",*t1,&n,&n,pivot,*green_plus,&n,&n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) t1[i][j]=green_plus[i][j];
		zgemm ("n", "n", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm ("c", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
#endif
		
		error = 0.0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				t1[i][j] = 0.5*(t1[i][j] + sigma_right[i][j]);
				a = real(t1[i][j]) - real(sigma_right[i][j]);
				b = imag(t1[i][j]) - imag(sigma_right[i][j]);
				error += sqrt (a * a +b * b);
				sigma_right[i][j] = t1[i][j];
			}
		error /= n * n;
	}
	
	
	
	/* left self-energy */  
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
          		sigma_left[i][j] = zero;
		sigma_left[i][i] = Complex (0.0, -1.0);
	}
	
	error = 1.0;
	for (iter = 0; iter < maxiter && error > maxerror; iter++) {

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
          			if (i!=j) 
						t1[i][j] = 0.0 - h0[i][j];
					else 
						t1[i][j] = e - h0[i][j];
					t1[i][j] = t1[i][j]- sigma_left[i][j];
          			
			}
#ifndef _AIX	
		zgetrf_ (&n, &n, *t1, &n, pivot, &x0);
		zgetri_ (&n, *t1, &n, pivot, *t0, &m, &x0);
		zgemm_ (&c123,&c123, &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm_ (&c123,&c123, &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
#else
		for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) green_plus[i][j]=zero;
			green_plus[i][i]=one; 
		}
		zgef(*t1,&n,&n,pivot);
		zgesm("n",*t1,&n,&n,pivot,*green_plus,&n,&n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) t1[i][j]=green_plus[i][j];
		/* t1 = h1^H * t1 * h1 */

		zgemm ("n", "c", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm ("n", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
		
		
#endif	
		error = 0.0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				t1[i][j] = 0.5*(t1[i][j] + sigma_left[i][j]);
				
				a = real(t1[i][j]) - real(sigma_left[i][j]);
				b = imag(t1[i][j]) - imag(sigma_left[i][j]);
				error += sqrt (a*a + b*b);
				sigma_left[i][j] = t1[i][j];
			}
		error /= n * n;
	}
	
	
	/* construct G+ */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
          		green_plus[i][j] = zero;
		}
		green_plus[i][i] = Complex (mu, 0.0);
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			green_plus[i][j]=green_plus[i][j]-h0[i][j];
		}

	/* self energies */
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
          	green_plus[i][j]=green_plus[i][j]-sigma_right[i][j];
			green_plus[i][j]=green_plus[i][j]-sigma_left[i][j];
          	
		}
#ifndef _AIX
	zgetrf_ (&n, &n, *green_plus, &n, pivot, &x0);
	if (x0)
		cerr<<"FATAL -- INTERNAL --\n"
			"Problems factorizing green_plus\n"
			"Stopping execution at "<< __FILE__<<"  "<<__LINE__;
	
	zgetri_ (&n, *green_plus, &n, pivot, *t0, &m, &x0);
      	if (x0)
		cerr<<"FATAL -- INTERNAL --\n"
			"Problems inverting green_plus\n"
			"Stopping execution at "<< __FILE__<<" "<< __LINE__;
#else
	for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) t1[i][j]=zero;
		t1[i][i]=one; 
	}
	zgef(*green_plus,&n,&n,pivot);
	zgesm("n",*green_plus,&n,&n,pivot,*t1,&n,&n);
	for (i=0;i<n;i++) for (j=0;j<n;j++) green_plus[i][j]=t1[i][j];
	
#endif	
		
		
	/* construct G+ - G- */
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++) {
			 error= 0.5 * (imag(green_plus[j][i])
						+ imag(green_plus[i][j]));
			green_plus[j][i] = Complex(error,0.5 * (-real(green_plus[j][i])
						+ real(green_plus[i][j])));
			
		}
	
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++) {
			green_plus[i][j] =  Complex(real(green_plus[j][i]),
			 -imag(green_plus[j][i]));
		}
	
	
	/* Kubo formula */
#ifndef _AIX
	zgemm_ (&c123,&c123, &n, &n, &n, &one, *x, &n, *green_plus, &n,
			&zero, *t0,  &n);
	zgemm_ (&c123,&c123, &n, &n, &n, &one, *t0, &n, *t0, &n, &zero,
			*green_plus, &n);
#else
	/* zgemm(*v,&n,"n",*green_plus,&n,"n",*t0,&n,&n,&n,&n,NULL,0);
	zgemm(*t0,&n,"n",*t0,&n,"n",*green_plus,&n,&n,&n,&n,NULL,0); */
	zgemm ("n", "n", &n, &n, &n, &one, *x, &n, *green_plus, &n,
			&zero, *t0,  &n);
	zgemm ("n", "n", &n, &n, &n, &one, *t0, &n, *t0, &n, &zero,
			*green_plus, &n);
#endif
	cond = 0.0;
      	for (i = 0; i < n; i++) {
       		cond +=  real(green_plus[i][i]);
	}
	return (2.0 *cond);
}
#endif // model_kondo_2n

#ifdef MODEL_KONDO_1N
double conductance3d(MyMatrix<MatType> const & matrix,double mu, int n,int dim,int l1,int maxiter,double maxerror)
{
	static	int	firstcall = 1;
	static int volume, *pivot;
	static	Complex	**h0, **h1, **t0, **t1, **x;
	static	Complex	**sigma_right, **sigma_left, **green_plus;
	/* tpem_sparse	*matrix; */
	int	 m, i, j, x0, x1, y0, y1, z0, z1, iter;
	unsigned int k;
	double	eta = 0.0, a,b, error, cond;
	Complex	one = Complex (1.0, 0.0), zero = Complex (0.0, 0.0), e;
	
	if (n == 0) {
		if (pivot) delete [] pivot, pivot = NULL;
		matrix_free (h0), h0 = NULL;
		matrix_free (h1), h1 = NULL;
		matrix_free (t0), t0 = NULL;
		matrix_free (t1), t1 = NULL;
		matrix_free (x), x = NULL;
		matrix_free (sigma_right), sigma_right = NULL;
		matrix_free (sigma_left),  sigma_left  = NULL;
		matrix_free (green_plus),  green_plus  = NULL;
		firstcall = 1;
		return -1;
	}
	
	/* matrix = (tpem_sparse *) this->data; */
	m = n * n;
	if (firstcall) {
		firstcall = 0;
		/* dimension = *(int *) parser_get ("dimension");
		volume    = *(int *) parser_get ("volume"); */
		pivot = new int [n]; 
		volume=n;
		
		h0 = matrix_alloc (n, n);
		h1 = matrix_alloc (n, n);
		t0 = matrix_alloc (n, n);
		t1 = matrix_alloc (n, n);
		x  = matrix_alloc (n, n);
		sigma_right = matrix_alloc (n, n);
		sigma_left  = matrix_alloc (n, n);
		green_plus  = matrix_alloc (n, n);
	}
	
	/* intralayer (resp. interlayer) hamiltonian h0 (resp. h1) */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			x[i][j] = zero;
			h0[i][j]=h1[i][j]=0;
		}
	}
	
	
	
	switch (dim) {
	case 3:
		for (i = 0; i < l1; i++) {
		for (j = 0; j < l1; j++) {
			for (k=0;k<l1;k++) {
				x0=k*l1*l1+j*l1+i;
				x[x0][x0]=k;
			}
		}
	}
	
	for (i = 0; i < l1; i++) {
	for (j = 0; j < l1; j++) {
	for (k = 0; k < l1; k++) {
		x1 = k*l1*l1 + j*l1 +i;
		for (x0=0;x0<l1;x0++) {
			for (y0=0;y0<l1;y0++) {
				for (z0=0;z0<l1;z0++) {
					y1=z0*l1*l1+y0*l1+x0;
					if (z0==0 && k==l1-1) {
						h1[x1][y1]=matrix(x1,y1);
					}
					else if (k!=0 || z0!=l1-1) {
						h0[x1][y1]=matrix(x1,y1);
					}
				}
			}
		}
	}}}
	break;
	case 2:
		for (i = 0; i < l1; i++) {
		for (j = 0; j < l1; j++) {
			x0=j*l1+i;
			x[x0][x0]=j;
		}
		}
	
	for (i = 0; i < l1; i++) {
	for (j = 0; j < l1; j++) {
		x1 =j*l1 +i;
		for (x0=0;x0<l1;x0++) {
			for (y0=0;y0<l1;y0++) {
					y1=y0*l1+x0;
					if (y0==0 && j==l1-1) {
						h1[x1][y1]=matrix(x1,y1);
					}
					else if (j!=0 || y0!=l1-1) {
						h0[x1][y1]=matrix(x1,y1);			
					}
				}
			}
	}}
	break;
	case 1:
	cerr<<"Cannot calculate conductance in 1d\n";
	return -1;
	break;
}
	
#ifndef _AIX
	/* velocity operator */
	zgemm_ ("n", "n", &n, &n, &n, &one, *h0, &n, *x, &n, &zero, *t0, &n);
    zgemm_ ("n", "n", &n, &n, &n, &one, *x, &n, *h0, &n, &zero, *t1, &n);
#else
	/* zgemm(*h0,&n,"n",*x,&n,"c",*t0,&n,&n,&n,&n,NULL,0);
	zgemm(*x,&n,"n",*h0,&n,"n",*t1,&n,&n,&n,&n,NULL,0); */
	zgemm ("n", "n", &n, &n, &n, &one, *h0, &n, *x, &n, &zero, *t0, &n);
    zgemm ("n", "n", &n, &n, &n, &one, *x, &n, *h0, &n, &zero, *t1, &n);
#endif
	
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
        		x[i][j] = Complex((double)(imag(t0[i][j] - t1[i][j])) /
				(l1 - 1.0),
		(double)(real(t1[i][j] - t0[i][j])) /
				(l1 - 1.0));
		}
	
	e = Complex (mu, eta);

	/* right self-energy */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
          		sigma_right[i][j] = zero;
		sigma_right[i][i] = Complex (0.0, -1.0);
	}
	
	error = 1.0;
	for (iter = 0; iter < maxiter && error > maxerror; iter++) {

		for (i = 0; i < n; i++) 
			for (j = 0; j < n; j++) {
				if (i!=j)
          			t1[i][j] =  0.0- h0[i][j]- sigma_right[i][j];
				else
					t1[i][j] = e - h0[i][j]- sigma_right[i][j];
          			
			}
#ifndef _AIX	
		zgetrf_ (&n, &n, *t1, &n, pivot, &x0);
		zgetri_ (&n, *t1, &n, pivot, *t0, &m, &x0);
		/* t1 = h1 * t1 * h1^H */
		zgemm_ ("n", "n", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm_ ("c", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
	
#else
		for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) green_plus[i][j]=zero;
			green_plus[i][i]=one; 
		}
		zgef(*t1,&n,&n,pivot);
		zgesm("n",*t1,&n,&n,pivot,*green_plus,&n,&n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) t1[i][j]=green_plus[i][j];
		zgemm ("n", "n", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm ("c", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
#endif
		
		error = 0.0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				t1[i][j] = 0.5*(t1[i][j] + sigma_right[i][j]);
				a = real(t1[i][j]) - real(sigma_right[i][j]);
				b = imag(t1[i][j]) - imag(sigma_right[i][j]);
				error += sqrt (a * a +b * b);
				sigma_right[i][j] = t1[i][j];
			}
		error /= n * n;
	}
	
	
	/* left self-energy */  
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
          		sigma_left[i][j] = zero;
		sigma_left[i][i] = Complex (0.0, -1.0);
	}
	
	error = 1.0;
	for (iter = 0; iter < maxiter && error > maxerror; iter++) {

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
          			if (i!=j) 
						t1[i][j] = 0.0 - h0[i][j];
					else 
						t1[i][j] = e - h0[i][j];
					t1[i][j] = t1[i][j]- sigma_left[i][j];
          			
			}
#ifndef _AIX	
		zgetrf_ (&n, &n, *t1, &n, pivot, &x0);
		zgetri_ (&n, *t1, &n, pivot, *t0, &m, &x0);
		zgemm_ ("n", "c", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm_ ("n", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
#else
		for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) green_plus[i][j]=zero;
			green_plus[i][i]=one; 
		}
		zgef(*t1,&n,&n,pivot);
		zgesm("n",*t1,&n,&n,pivot,*green_plus,&n,&n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) t1[i][j]=green_plus[i][j];
		/* t1 = h1^H * t1 * h1 */

		zgemm ("n", "c", &n, &n, &n, &one, *t1, &n, *h1, &n,
				&zero, *t0, &n);
	    zgemm ("n", "n", &n, &n, &n, &one, *h1, &n, *t0, &n,
				&zero, *t1, &n);
		
		
#endif	
		error = 0.0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				t1[i][j] = 0.5*(t1[i][j] + sigma_left[i][j]);
				
				a = real(t1[i][j]) - real(sigma_left[i][j]);
				b = imag(t1[i][j]) - imag(sigma_left[i][j]);
				error += sqrt (a*a + b*b);
				sigma_left[i][j] = t1[i][j];
			}
		error /= n * n;
	}
	
	if (error>maxerror) {
		//cerr<<__FILE__<<" "<<__LINE__<<" in vCond error is "<<error;
		//cerr<<" but maxerror="<<maxerror<<endl;
		//exit(1);
		return -100.;
	}
	
	
	/* construct G+ */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
          		green_plus[i][j] = zero;
		}
		green_plus[i][i] = Complex (mu, 0.0);
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			green_plus[i][j]=green_plus[i][j]-h0[i][j];
		}

	/* self energies */
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
          	green_plus[i][j]=green_plus[i][j]-sigma_right[i][j];
			green_plus[i][j]=green_plus[i][j]-sigma_left[i][j];
          	
		}
#ifndef _AIX
	zgetrf_ (&n, &n, *green_plus, &n, pivot, &x0);
	if (x0)
		cerr<<"FATAL -- INTERNAL --\n"
			"Problems factorizing green_plus\n"
			"Stopping execution at "<< __FILE__<<"  "<<__LINE__;
	
	zgetri_ (&n, *green_plus, &n, pivot, *t0, &m, &x0);
      	if (x0)
		cerr<<"FATAL -- INTERNAL --\n"
			"Problems inverting green_plus\n"
			"Stopping execution at "<< __FILE__<<" "<< __LINE__;
#else
	for (i=0;i<n;i++) { 
		for (j=0;j<n;j++) t1[i][j]=zero;
		t1[i][i]=one; 
	}
	zgef(*green_plus,&n,&n,pivot);
	zgesm("n",*green_plus,&n,&n,pivot,*t1,&n,&n);
	for (i=0;i<n;i++) for (j=0;j<n;j++) green_plus[i][j]=t1[i][j];
	
#endif	
		
		
	/* construct G+ - G- */
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++) {
			 error= 0.5 * (imag(green_plus[j][i])
						+ imag(green_plus[i][j]));
			green_plus[j][i] = Complex(error,0.5 * (-real(green_plus[j][i])
						+ real(green_plus[i][j])));
			
		}
	
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++) {
			green_plus[i][j] =  Complex(real(green_plus[j][i]),
			 -imag(green_plus[j][i]));
		}
	
	
	/* Kubo formula */
#ifndef _AIX
	zgemm_ ("n", "n", &n, &n, &n, &one, *x, &n, *green_plus, &n,
			&zero, *t0,  &n);
	zgemm_ ("n", "n", &n, &n, &n, &one, *t0, &n, *t0, &n, &zero,
			*green_plus, &n);
#else
	/* zgemm(*v,&n,"n",*green_plus,&n,"n",*t0,&n,&n,&n,&n,NULL,0);
	zgemm(*t0,&n,"n",*t0,&n,"n",*green_plus,&n,&n,&n,&n,NULL,0); */
	zgemm ("n", "n", &n, &n, &n, &one, *x, &n, *green_plus, &n,
			&zero, *t0,  &n);
	zgemm ("n", "n", &n, &n, &n, &one, *t0, &n, *t0, &n, &zero,
			*green_plus, &n);
#endif
	cond = 0.0;
      	for (i = 0; i < n; i++) {
       		cond +=  real(green_plus[i][i]);
	}
	return (2.0 *cond);
}


			
#endif



