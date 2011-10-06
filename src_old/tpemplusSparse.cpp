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


tpem_sparse *new_tpem_sparse (unsigned int rank, unsigned int size)
{
	tpem_sparse *t = new tpem_sparse;
	

	t->rank = rank;
	t->rowptr = new  unsigned int[rank + 1];
	t->colind = new  unsigned int[size];
	t->values = new  tpem_t[size];
	t->rowptr[rank] = size;
	return t;
}


void tpem_sparse_copy (tpem_sparse *t, tpem_sparse *other)
{
	unsigned int i, n;

	n = other->rank = t->rank;
	for (i = 0; i < n; i++)
		other->rowptr[i] = t->rowptr[i];
	n = other->rowptr[n] = t->rowptr[n];
	for (i = 0; i < n; i++) {
		other->colind[i] = t->colind[i];
		other->values[i] = t->values[i];
	}
}


tpem_t *tpem_sparse_element (tpem_sparse *t, unsigned int row, unsigned int col)
{
	unsigned int i;

	for (i = t->rowptr[row]; i < t->rowptr[row + 1]; i++)
		if (t->colind[i] == col)
			return t->values + i;
	cerr<<"tpem_sparse_element: no element at ["<<row<<"]["<<col<<"]\n";
	cerr<<"tpem_sparse_element: debug: rowprtr[row]="<<t->rowptr[row];
	cerr<<",rowptr[row+1]="<<t->rowptr[row+1]<<endl;
	return NULL;
}

void tpem_sparse_free (tpem_sparse *t)
{
	delete [] t;
	t=0;
}

/*
tpem_t **new_tpem_matrix (unsigned int row, unsigned int col)
{
	unsigned int i;
	tpem_t **t = new tpem_t *[row];

	*t = new tpem_t[row * col];
	for (i = 1; i < row; i++)
		t[i] = t[i - 1] + col;
	return t;
}
*/
