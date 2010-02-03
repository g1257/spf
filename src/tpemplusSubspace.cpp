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


void tpem_subspace_free (tpem_subspace *t)
{
	if (t && t->stack) delete [] (t->stack), t->stack = NULL;
	if (t && t->flags) delete [] (t->flags), t->flags = NULL;
	if (t) delete [] (t), t = NULL;
}


void tpem_subspace_reset (tpem_subspace *t)
{
	unsigned int i;

	for (i = 0; i < t->size; i++)
		t->flags[i] = 0;
	t->top = t->stack;
}


tpem_subspace *new_tpem_subspace (unsigned int size)
{
	tpem_subspace *t = new tpem_subspace;

	t->size = size;
	t->flags =new int [size];
	t->stack =new unsigned int [size]; 
	tpem_subspace_reset (t);
	return t;
}

void tpem_subspace_fill(tpem_subspace *t)
{
	for (size_t state=0;state<t->size;state++) tpem_subspace_push(t,state);
}

void tpem_subspace_push (tpem_subspace *t, unsigned int state)
{
	if (t->flags[state] == 0) {
		t->flags[state] = 1;
		*(t->top)++ = state;
	}
}
