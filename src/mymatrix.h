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




#ifndef MYMATRIX_H_
#define MYMATRIX_H_
#include <iostream>
#include <vector>
#include "basic_types.h"


template<class T> class MyMatrix;  // pre-declare the template class itself
template<class T> std::ostream& operator<< (std::ostream& o, MyMatrix<T> const & x); 

using std::allocator;



/** \brief Implements a full n x n matrix. Used only by DIAG algorithm */
template<class T>
class MyMatrix {
	std::vector<T> data;
	int n;
	int total;
	void copy(std::vector<T> &array) const; //auxiliary member function
public:
		MyMatrix() { n=0; }
                ~MyMatrix();   
                void init(int,T=0);

		T operator()(int i,int j) const { 
			return data[i+j*n];
		}
		T &operator()(int i,int j) {
			return data[i+j*n];
		}
		void set(int i,int j,T value) {
			data[i+j*n]=value;
		}
		MyMatrix &operator= (const MyMatrix& myObject); // assignment operator
		MyMatrix(MyMatrix const &myObject); // copy constructor
		int getRank() const { return n; }
		int isHermitian() const;
		friend std::ostream &operator<< <>(std::ostream &s,MyMatrix<T> const &myObject);
};

int diag(MyMatrix<MatType> &myMatrix,std::vector<double> &e,char jobz='n');

#endif
