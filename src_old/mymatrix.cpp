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




#include "mymatrix.h"

#ifdef __linux
#ifndef __sgi
#include "diag_linux.h"
#endif
#endif

#ifdef _AIX
#include "diag_aix.h"
#endif

#ifdef _CRAY
#include "diag_cray.h"
#endif

#ifdef __macosx
#include "diag_macosx.h"
#endif

#ifdef __sgi
#include "diag_sgi.h"
#endif

template<class T>
int MyMatrix<T>::isHermitian() const
{
	int i,j;
	T v;
	double tmp;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			v=data[i+j*n]-conj(data[j+i*n]);
			tmp = real(v)*real(v)+imag(v)*imag(v);
			if (tmp>1e-8) {
				std::cerr<<"Not hermitian when i="<<i<<" and j="<<j;
				std::cerr<<" "<<data[i+j*n]<<"!="<<data[j+i*n];
				std::cerr<<" value is "<<tmp<<std::endl;
				return 0;
			}
		}
	}
	return 1;
}

template<class T>
void MyMatrix<T>::init(int n1,T value)
{
        n = n1;
        total = n*n;
        data.insert(data.begin(),n*n,value);
}


template<class T>
MyMatrix<T>::~MyMatrix()
{
	data.clear();
}

template<class T>
MyMatrix<T> &MyMatrix<T>::operator=(const MyMatrix &myMatrix)
{
	if (this == &myMatrix) return *this;   // Gracefully handle self assignment[12.1]

	// Put the normal assignment duties here, which in this case would be:
	n = myMatrix.getRank();
	data.clear();
	data.insert(data.begin(),n*n,0);
	myMatrix.copy(data);
	return *this;
}

template<class T>
void MyMatrix<T>::copy(std::vector<T> &array) const { //auxiliary member function
	for (int i=0;i<total;i++) array[i]=data[i]; 
}


template<class T>
std::ostream &operator<<(std::ostream &s,MyMatrix<T> const &myMatrix)
{
	int i,j;
	for (i=0;i<myMatrix.n;i++) {
		for (j=0;j<myMatrix.n;j++) {
			s<<myMatrix(i,j)<<" ";
		}
		s<<std::endl;
	}
	return s;
}




// The line below is needed because it is not (directly) possible to
// separate a template class in .h and .cpp parts
template class MyMatrix<MatType>;
template
std::ostream &operator<< (std::ostream &s,MyMatrix<MatType> const &myMatrix);
