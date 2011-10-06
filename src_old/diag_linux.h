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




#ifndef NO_LAPACK
extern "C" void zheev_(char &,char &,int &,MatType [],int &,
                       double [],MatType [],int &,double [],int &);
/* extern "C" void zgemm_(char &,char &,int &,int &,int &,MatType &,
		MatType [],int &,MatType [],int &,MatType &,MatType [],int &);
*/

int diag(MyMatrix<MatType> &matrix,std::vector<double> &eig,char jobz)
{
	int info,lwork,i,j;
        char uu='U',jobz2;
        int matsize = matrix.getRank();
        //static int firstcall=1;        
        lwork = 3*matsize;
        static std::vector<MatType> work(lwork+10);
        static std::vector<double> rwork(3*matsize-1);
        
        /*work = new MatType[lwork+10];
        rwork = new double[3*matsize-1];
*/
        jobz2=jobz;
        zheev_(jobz2,uu,matsize,&(matrix.pointer()[0]),matsize,&(eig[0]),&(work[0]),lwork,&(rwork[0]),info);

        /*delete [] rwork;
        delete [] work;*/
        return info;
}
#else
extern "C" {
#include "diag.h"
}

int diag(MyMatrix<MatType> &matrix,std::vector<double> &eig,char jobz)
{
	ccomplex *a;
	unsigned int matsize = matrix.getRank();
	double *e;
	static int firstcall=1;
	unsigned int i,j;
	
	
	a = new ccomplex[matsize * matsize];
	e = new double[matsize];
	
	for (i=0;i<matsize;i++) {
                for (j=0;j<matsize;j++) {
			a[i+j*matsize].re = real(matrix(i,j));
			a[i+j*matsize].im = imag(matrix(i,j));
		}
	}
	diag(matsize,a,e);
	for (i=0;i<matsize;i++) eig[i]=e[i];
	delete [] e;
	delete [] a;
	return 1;
}	
		
#endif

		
