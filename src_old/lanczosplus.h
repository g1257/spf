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




#ifndef LANCZOSPLUS_H_
#define LANCZOSPLUS_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "basic_types.h"
#include "tpemplusTypes.h"

using std::vector;


class SparseMatrix {
    public:

 /* constructor, takes a tpem_sparse to create a SparseMatrix */  
  SparseMatrix(tpem_sparse const *what)
  {
    MatType c;
    int i;
    size = what->rank;
    
    
	for (i=0;i<size+1;i++) rowptr.push_back(what->rowptr[i]);
	
	for (i=0;i<rowptr[size];i++) colind.push_back(what->colind[i]);
	
	
    for (i=0; i<rowptr[size]; i++)
      {
#ifdef TPEM_COMPLEX
	c = what->values[i];
	values.push_back(c);
	
#else
	c = MatType(what->values[i], 0.0);
	values.push_back(c);
		
#endif
      }
  } // end of constructor

            ~SparseMatrix() { rowptr.clear(); colind.clear(); values.clear();  }
            void initRow(int n) {
                rowptr.insert(rowptr.begin(),n);
            }
            void initCol(int n) {
                colind.insert(colind.begin(),n);
            }
            void initValues(int n) {
                values.insert(values.begin(),n);
            }
            void setRow(int n,int v) {
                rowptr[n]=v;
            }
            void setCol(int n,int v) {
                colind[n]=v;
            }
            void setValues(int n,MatType &v) { 
                values[n]=v;
            }    
            
	    void debug() const {
	    	 std::cerr<<"debug rowptr[0]="<<rowptr[0]<<std::endl;
	}
		
            void sparse_mult (vector<MatType> &x, vector<MatType> const &y) const
            { 
                   /* performs x = x + A * y
                   * where x and y are vectors and A is a sparse matrix in
                   * row-compressed format */
                  unsigned int i;
		  int j;
		
              for (i = 0; i < y.size(); i++) {
	      	
                for (j = rowptr[i]; j < rowptr[i + 1]; j++) {
                  x[i] += values[j] * y[colind[j]];
		 
		 }
		}
		 

            }
            int getSize() const { return size; }
	    void setSize(int s) { size = s; }
                        
    private:
            std::vector<int> rowptr;
            std::vector<int> colind;
            std::vector<MatType> values;
	    int size;
            
}; // sparsematrix

// export this function only
int lanczos (int flag, SparseMatrix const &mat, int max_nstep, double &gsEnergy, std::vector<MatType> &z);

#endif

