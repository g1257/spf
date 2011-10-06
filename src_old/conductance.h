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




#ifndef MESCOND_H
#define MESCOND_H
#include "basic.h"
#include "dynvars.h"
#include "geometry.h"
#include "mymatrix.h"
#include "dynvars.h"
#include "aux.h"
#include "parameters.h"

typedef  std::complex<double> Complex;

#ifndef _AIX
        extern  "C" void        zgemm_  (char * , char *, int *, int *, int *,
                                Complex *, Complex *, int *, Complex *, int *,
                                Complex *, Complex *, int *);

        extern  "C" void        zgetrf_ (int *, int *, Complex *, int *, int *, int *);
        extern  "C" void        zgetri_ (int *, Complex *, int *, int *, Complex *,
                                int *, int *);
#else
	extern "C" void zgemm (char *,char *,int *,int *, int *,Complex *,Complex *,int *,
	Complex *,int *,Complex *,Complex *,int *);
	extern "C" void zgef (Complex *,int *,int *,int *);
	extern "C" void zgesm (char *,Complex *,int *,int *,int *,Complex *,int *,int *);
#endif


#ifndef _AIX
        extern  "C" void        zgemm_  (char * , char *, int *, int *, int *,
                                Complex *, Complex *, int *, Complex *, int *,
                                Complex *, Complex *, int *);

        extern  "C" void        zgetrf_ (int *, int *, Complex *, int *, int *, int *);
        extern  "C" void        zgetri_ (int *, Complex *, int *, int *, Complex *,
                                int *, int *);
#endif


extern void setupHamiltonian(MyMatrix<MatType> & matrix,Geometry const &geometry, DynVars const &dynVars, 
        Parameters const &ether,int type);
extern Complex **matrix_alloc (unsigned int row, unsigned int col);
extern void matrix_free (Complex **this1a);	
double conductance3d(MyMatrix<MatType> const & matrix,double mu, int n,
int dim,int l1,int maxiter,double maxerror);
	
#endif
