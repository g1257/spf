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




/***************************************************************************************************
* TpemPlus Library by G. A.
* 
****************************************************************************************************/
#include "tpemplus.h"

extern "C" {
#include "diag.h"
#ifdef  TPEM_COMPLEX
void    zheev_ (char *, char *, int *, tpem_t *, int *,
                        double *, tpem_t *, int *, double *, int *);
#else
void    dsyev_ (char *, char *, int *, tpem_t *, int *,
                        double *, tpem_t *, int *, int *);
#endif
}

double my_random ()
{
	static int next = 1;
	
	next = 16807 * (next % 127773) - 2836 * (next / 127773);
	if (next <= 0) next += 2147483647;
	return ((double) next) / 2147483647.0;
}


tpem_sparse *new_tpem_sparse_random (unsigned int rank, double range)
{
	tpem_sparse *t = new_tpem_sparse (rank, 3 * rank);
	double scale = 2.0 + fabs (range), epsilon;
	tpem_t	hopping = -1.0 / scale;
	unsigned int	i, j;

	range *= 2.0 / scale;
	for (i = j = 0; i < rank; i++) {
		t->rowptr[i] = j;
		t->colind[j] = i;
		epsilon = range * (my_random () - 0.5);
		t->values[j] = epsilon;
		j++;
		t->colind[j] = (i + 1) % rank;
		t->values[j] = hopping;
		j++;
		t->colind[j] = (i + rank - 1) % rank;
		t->values[j] = hopping;
		j++;
	}
	return t;
}


double tpem_apply (tpem_sparse *matrix, tpem_func_ptr funcptr,TpemOptions const &tpemOptions)
{
	int cutoff=tpemOptions.cutoff;
	vector<double>	coeffs(cutoff),moment(cutoff);
	double	ret;
	
	tpem_calculate_coeffs (coeffs, funcptr,tpemOptions);
	tpem_calculate_moment(matrix, moment,tpemOptions);
	ret = tpem_expansion (moment, coeffs);

	return ret;
}

#ifndef NO_LAPACK
/* Calculates an observable defined by the function funcptr using exact diagonalization
    for the model given by the Hamiltonian matrix *matrix */

static double diag_apply (tpem_sparse *matrix, double (*funcptr)(double))
{
	char	jobz	= 'n';
	char	uplo	= 'u';
	int	n	= matrix->rank, info, i;
	unsigned int	k;
#ifdef	TPEM_COMPLEX
	int	lwork 	= 2 * n - 1;
	double	*rwork	= new double[(3 * n - 2)];
#else
	int	lwork	= 3 * n - 1;
#endif
	tpem_t	*a	= new tpem_t[n * n];
	double	*w	= new double[n];
	double ret;
	tpem_t	*work	= new tpem_t[lwork]; 

	for (i = 0; i < n * n; i++)
		a[i] = 0.0;
	for (i = 0; i < n; i++)
		for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++)
			a[i * n + matrix->colind[k]] = matrix->values[k];
#ifdef	TPEM_COMPLEX
	zheev_ (&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
	delete [] rwork;
#else
	dsyev_ (&jobz, &uplo, &n, a, &n, w, work, &lwork, &info);
#endif
	ret = 0.0;
	for (i = 0; i < n; i++) {
		ret += funcptr (w[i]);
	}	
	delete [] a;
	delete [] w;
	delete [] work;
	return ret;
}

#else
static double diag_apply (tpem_sparse *matrix, double (*funcptr)(double))
{

        int     n       = matrix->rank, i;
        unsigned int  k;
        ccomplex *a;
        double ret;
        double *w;

        a = new ccomplex [n*n];
        w = new double[n];
        for (i = 0; i < n * n; i++)
                a[i].re = a[i].im = 0.0;
        for (i = 0; i < n; i++) {
                for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++) {
                        a[i * n + matrix->colind[k]].re = real(matrix->values[k]);
                        a[i * n + matrix->colind[k]].im = imag(matrix->values[k]);
                }
        }
        diag(n, a, w);

        ret = 0.0;
        for (i = 0; i < n; i++) {
                ret += funcptr (w[i]);
        }
        delete [] a;
        delete [] w;

        return ret;
}
#endif

double tpem_apply_diff (tpem_sparse *matrix0, tpem_sparse *matrix1,
		tpem_func_ptr funcptr, TpemOptions const &tpemOptions)
{
	int cutoff=tpemOptions.cutoff;
	vector<double> coeffs(cutoff);
	vector<double> moment(cutoff);
	vector<unsigned int> support(2);
	double  ret;

	support[0] = 0;
	support[1] = matrix0->rank / 2 - 1;
	tpem_calculate_coeffs (coeffs, funcptr,tpemOptions);

	tpem_calculate_moment_diff (matrix0,matrix1,  moment,support,tpemOptions);
	ret = tpem_expansion (moment, coeffs);
	
	return ret;
}



double orig_func1(double x)
{
	return 5.0 * x *  (1.0 - tanh (10.0 * x));
}

double orig_func2(double x)
{
	return 0.5 * (1.0 - tanh (10.0 * x));
}

double orig_func3(double x)
{
	return log (1.0 + exp (-10.0 * x));
}



using namespace std;

int main (int argc,char *argv[])
{
	tpem_sparse *matrix0 = new_tpem_sparse_random (400, 10.0);
	tpem_sparse *matrix1 = new_tpem_sparse_random (400, 10.0);
	double	tpem0, tpem1, tpemd, naive0, naive1;
	unsigned int	cutoff;
	int rank;
	TpemOptions tpemOptions;
	
	tpemOptions.mpi_nop2=1;
        tpem_init(argc,argv,tpemOptions);
	rank = tpemOptions.rank;
	tpemOptions.epsProd=1e-5;
	tpemOptions.epsTrace=1e-7;
	tpemOptions.tpemType="tpem";
	tpemOptions.coeffs = "accurate";
	
	if (rank==0) {
	cout<<"********************************************************\n";
	cout<<"****** TESTING TRUNCATED POLYNOMIAL EXPANSION **********\n";
	cout<<"********************************************************\n";
	cout<<"\n";
	cout<<"\n";
	cout<<"this testing program calculates model properties in two ways:\n";
	cout<<"(i)  Using standard diagonalization\n";
	cout<<"(ii) Using the truncated polynomial expansion method\n";
	cout<<"\n";
	cout<<"All tests are done for a nearest neighbour interaction with\n";
	cout<<"random (diagonal) potentials.\n";
	cout<<"\n";
	cout<<"\n";
	}
	
	tpem_sparse_copy (matrix0, matrix1);
	matrix1->values[matrix1->rowptr[0]] = 2.4 * (my_random () - 0.5);
	matrix1->values[matrix1->rowptr[matrix1->rank / 2 - 1]] = 2.4 * (my_random () - 0.5);

if (rank==0) {	
	cout<<"-------------------------------------------------------------\n";
	cout<<"TEST 1: MEAN VALUE FOR THE FUNCTION:                         \n";
	cout<<"        E(x) = 5.0 * x * (1.0 - tanh (10.0 * x))\n";
	cout<<"\n";
}
	naive0 = diag_apply (matrix0, orig_func1);
if (rank==0) {
	cout<<"** Using diagonalization <E>="<<naive0<<"\n";
	cout<<"** Using TPEM <E>=(cutoff--> infinity) "
			"lim<E_cutoff> where <E_cutoff> is\n";
	cout<<"cutoff\t<E_cutoff>\tError(compared to diag.)\n";
}
	for (cutoff = 20; cutoff <= 40; cutoff++) {
		tpemOptions.cutoff=cutoff;
		tpem0 = tpem_apply (matrix0, orig_func1, tpemOptions);
if (rank==0) {
		cout<<cutoff<<"\t"<<tpem0<<"\t"<<naive0<<"\t";
		cout<<(100.0 * fabs (1.0 - tpem0 / naive0))<<endl;
		cout.flush();
}
	}

if (rank==0) {	
	cout<<"-------------------------------------------------------------\n";
	cout<<"TEST 2: MEAN VALUE FOR THE FUNCTION:                         \n";
	cout<<"        N(x) =  0.5 * (1.0 - tanh (10.0 * x))";
	cout<<"\n";
}
	naive0 = diag_apply (matrix0, orig_func2);
if (rank==0) {
	cout<<"** Using diagonalization <N>="<<naive0<<endl;
	cout<<"** Using TPEM <N>=(cutoff--> infinity) "
			"lim<N_cutoff> where <N_cutoff> is \n";
	cout<<"cutoff\t<N_cutoff>\tError(compared to diag.)\n";
}
	for (cutoff = 20; cutoff <= 40; cutoff++) {
		tpemOptions.cutoff=cutoff;
		tpem0 = tpem_apply (matrix0, orig_func2, tpemOptions);
if (rank==0) {
		cout<<cutoff<<"\t"<<tpem0<<"\t";
		cout<<(100.0 * fabs (1.0 - tpem0 / naive0))<<endl;
		cout.flush();
}
	}

if (rank==0) {	
	cout<<"-------------------------------------------------------------\n";
	cout<<"TEST 3: MEAN VALUE AND DIFFERENCE FOR THE FUNCTION:          \n";
	cout<<"        S(x) =  log (1.0 + exp (-20.0 * x))\n";
	cout<<"\n";
}
	naive0 = diag_apply (matrix0, orig_func3);
	naive1 = diag_apply(matrix1, orig_func3);
if (rank==0) {
	cout<<"** Using diagonalization <S[matrix0]>= "<<naive0<<endl;
	cout<<"** Using diagonalization <S[matrix1]>= "<<naive1<<endl;
	cout<<"** Using diagonalization <S[matrix1]>-<S[matrix0]>="<<(naive0 - naive1)<<endl;
	cout<<"** Using TPEM <S>=(cutoff--> infinity) lim<S_cutoff>\n";
	cout<<"cutoff\tDelta_S_cutoff\tS_cutoff[diff]\tError (to diag.)\n";
}
	for (cutoff = 10; cutoff <= 40; cutoff++) {
		tpemOptions.cutoff=cutoff;
		tpem0 = tpem_apply (matrix0, orig_func3, tpemOptions);
		tpem1 = tpem_apply (matrix1, orig_func3, tpemOptions);
		tpemd =0.0;
		tpemd = tpem_apply_diff (matrix0, matrix1, orig_func3,tpemOptions);
if (rank==0) {
		cout<<cutoff<<"\t"<<(tpem0-tpem1)<<"\t"<<tpemd<<"\t";
		cout<<(100.0 * fabs (1.0 - tpemd / (naive0-naive1)))<<endl;
		cout.flush();
}
	}

   	
	tpem_sparse_free (matrix0);
	tpem_sparse_free (matrix1);
	
	tpem_finalize();
	return 0;
}
