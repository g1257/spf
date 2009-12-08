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




extern "C" {
#include <vecLib/vecLib.h>
#include <vecLib/clapack.h>
}

int diag(MyMatrix<MatType> &matrix,std::vector<double> &eig,char jobz)
{
	int info,lwork,i,j;
	__CLPK_doublecomplex *work,*a;
	double *rwork,*ee;
	char uu='U',jobz2;
	int matsize = matrix.getRank();
	static int firstcall=1;
	MatType tmp;

	lwork = 3*matsize;
	
	work = new __CLPK_doublecomplex[lwork+10];
	rwork = new double[3*matsize-1];
	a  = new __CLPK_doublecomplex[matsize*matsize];
	ee = new double[matsize];
	
	for (i=0;i<matsize;i++) {
		ee[i]=eig[i];
		for (j=0;j<matsize;j++) {
			a[j+i*matsize].r=real(matrix(j,i));
			a[j+i*matsize].i=imag(matrix(j,i));
			//cout<<"matrix["<<j<<"]["<<i<<"]="<<matrix[j][i]<<endl;
		}
	}
	//for (i=0;i<matsize*matsize;i++) cout<<"linindex="<<i<<" "<<a[i];
	
	//cout<<"About to do diagonalization\n";
	jobz2=jobz;
	zheev_(&jobz2,&uu,&matsize,&(a[0]),&matsize,ee,work,&lwork,rwork,&info);
	//cout<<"Diagonalization done!!\n";
	for (i=0;i<matsize;i++) eig[i]=ee[i];
	
	if (jobz=='V' || jobz=='v') {
		for (i=0;i<matsize;i++) {
			for (j=0;j<matsize;j++) {
				tmp=MatType(a[i+j*matsize].r,a[i+j*matsize].i);
				matrix(i,j)=tmp;
			}
		}
	}
	delete [] rwork;
	delete [] work;
	delete [] a;
	delete [] ee;
	return info;

}
