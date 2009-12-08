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




#ifndef AUX_H_
#define AUX_H_
#include "tpemplusTypes.h"
#include "histogram.h"
#include "mymatrix.h"
/** \brief Auxiliary structures and parameters that change frequently. 

 * For parameters that don't usually change use the class Parameters.<br>
 * For dynamical variables use the class DynVars. This class has an
 * assignment operator to deal with the pointer data members, so
 * a thing like aux1 = aux2 is permitted. */

struct Aux {
	Histogram Nw[6];
	Histogram *Arw, *ArwC,*Ldos;
	Histogram Sigma;
	std::vector<double> eigM,lcd,clasCor,cco,oco;
	MyMatrix<MatType> matrix;
	std::vector<double> eigOneBand,eigAllBands;
	tpem_sparse **sparseMatrix;
	tpem_sparse **sparseTmp;
	std::vector<double> avMoments, curMoments;
	std::vector<MatType> offdMoments;
	std::vector<MatType> opticalMoments;
	// parameters that change
	double varMu, varTpem_a,varTpem_b,atmp,btmp;
	double varEnergy1, varEnergy2;
	vector<double> nac;
	vector<double> bcsCorxx,bcsCorxy;
	vector<double> orbitalAngles;
	
	/* Aux (Aux const &myObject) // copy constructor
	{
				
		mymatrix = myObject.mymatrix;
		eigOneBand = myObject.eigOneBand;
		eigAllBands = myObject.eigAllBands;
		
	} */
			
	Aux &operator= (const Aux& myObject)
	{
        	if (this == &myObject) return *this;   // Gracefully handle self assignment[12.1]

        	// Put the normal assignment duties here, which in this case would be:
        	matrix  = myObject.matrix;
		eigOneBand = myObject.eigOneBand;
		eigAllBands = myObject.eigAllBands;
		varMu = myObject.varMu;
		varTpem_a = myObject.varTpem_a;
		varTpem_b = myObject.varTpem_b;
		atmp = myObject.atmp;
		btmp = myObject.btmp;
		varEnergy1 = myObject.varEnergy1;
		varEnergy2 = myObject.varEnergy2;	
        	return *this;
	} 
};


#endif
