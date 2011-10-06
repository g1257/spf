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




#define zhpev esvzhpev
typedef union { struct { double _re, _im; } _data; double _align; } dcmplx;

extern "C" {
int esvzhpev(int, dcmplx *, double *, void *, int, int, double *, int );
}

int diag(MyMatrix<MatType> &matrix,std::vector<double> &eig,char jobz)
{
// Ok, AIX uses ESSL which is somewhat different from dxml or the 
// netlib version of lapack. The main difference is that the 
// matrix must be passed lower-packet for eigenvals only and
// uper-packet for eigenvects only (aix==aches)
	
	int i,j,k;	
	static MatType *work;
	static double *rwork;
	static int firstcall=1;
	int matsize = matrix.getRank();
	
	if (firstcall) {
		firstcall=0;
		work=new MatType[matsize*matsize];
		rwork = new double[matsize*4+1];
	}
	
	char ul = 'l';
	int lwork = 2*matsize-1;
	int info = 100;
	// Declare and allocate ap
	dcmplx *ap=new dcmplx[matsize*(matsize+1)/2];
	//MatType *ap = new MatType[matsize*(matsize+1)/2];

	// First pack the matrix into "ap"
	if (tolower(jobz)=='n') {
		k=0;
		for (j=0;j<matsize;j++) {
			for (i=j;i<matsize;i++) {
				ap[k]._data._re=real(matrix(i,j));
				ap[k++]._data._im=imag(matrix(i,j));
				//ap[k++]=matrix[i][j];
			}
		}
	}
	else {
		k=0;
		for (j=0;j<matsize;j++) {
			for (i=j;i<matsize;i++) {
				ap[k]._data._re=real(matrix(i,j));
				ap[k++]._data._im=imag(matrix(i,j));
				//ap[k++]=matrix[i][j];
			}
		}
	}
	
	// Now create space to store the eigenvectors
	
	
	if (tolower(jobz)=='n') i=0;
	else i=1;
	j=4*matsize;
	zhpev(i,ap,&(eig[0]),work,matsize,matsize,rwork,j);
	
	// Copy work into the original matrix mat 
	if (jobz=='V' || jobz=='v') {
		for (j=0;j<matsize;j++) {
			for (i=0;i<matsize;i++) {
				matrix.set(i,j,work[i+j*matsize]);
				// mat(i,j).im=work[i+j*dim].re;
			}
		}
	}
	
	// Deallocate everything
	delete [] ap;
	return info;
}
