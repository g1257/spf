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




/* An Histogram class
 *
 * by G. Alvarez (Feb, 2001)
 *
 * Constructor: Histogram sampleObject(minimumEnergy,maximumEnergy,deltaEnergy);
 *
 * Add a value: sampleObject.add(energyValue);
 *
 * Write histogram to strOut: sampleObject.print(strOut);
 *
 * See sample program for further details        
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <iostream>
#include <sstream>
#include <cstdlib>
		 		 
using std::ostream;
using std::cerr;
using std::endl;
		 

/*! \brief A class to deal with histograms

 *  - Constructor: Histogram sampleObject(minimumEnergy,maximumEnergy,deltaEnergy);
 *  - Add a value: sampleObject.add(energyValue);
 *  - Write histogram to strOut: sampleObject.print(strOut);
*/
class Histogram {
public:
	/** \brief Constructor
	
	 * \param minE_: Low energy bound
	 * \param maxE_: High energy bound
	 * \parm steps_: Number of histogram divisions.
	 */
	Histogram(double minE_,double maxE_,int 
		steps_):minE(minE_),maxE(maxE_),steps(steps_) 
	{
		// Some checking
    	 
		if (steps<=0 || minE>=maxE) {
			cerr<<"Histogram: steps: "<<steps<<" minE= "<<minE<<" maxE= "<<maxE<<endl;
			steps=0; // so we remember NOT to deallocate memory, see the destructor
			cerr<<"Histogram: FATAL AT: "<<__FILE__<<":"<<__LINE__<<endl;
			exit(1);
	         }
             
		// Memory allocation
		histE = new double[steps+1];
		histDE = new double[steps+1];
	         
		deltaE = (double)(maxE-minE)/steps;
			 
		for (int n=0;n<steps+1;n++) {
			histE[n]=minE+n*deltaE;
			histDE[n]=0;
		}    
	}      
	Histogram() { steps=0;	}
      
     
	~Histogram();
	/** Returns the number of steps */
	int size() const { return steps; }
	/** Returns the x-coordinate (center of interval) of the histogram at step i */  
	double coorX(int i) const  {return histE[i]; }
	/** Returns the y-coordinate (height) of the histogram at step i */
	double coorY(int i) const { return histDE[i]; }
	int dummy1() const { return 0; }
        
	/** Same as constructor. */
	void init(double minE_,double maxE_,int steps_,const double& initVal=0);
	
	/** Adds one to the step corresponding to energy Evalue */
	int add(double EValue); 
	/** Adds Evalue to step corresponding to energy xValue */
	int add(double xValue,double EValue);
	void dump(ostream & strOut) const { print(strOut); }
	/** Prints the histogram */
	void print(ostream & strOut) const {
		for (int n=0;n<steps;n++) strOut<<histE[n]<<" "<<histDE[n]<<endl;
	}

	int getSubInterval(const double& x) const;

	int multiply(double xValue,double EValue);

	/** Returns 1 if the histogram has been initalized with more than 0 steps or zero otherwise */
	int enabled() const {
		if (steps<=0) return 0;
		return 1;
	}	  
	  
	/** Returns the sum of all histogram values from energy e1 to energy e2 */
	double integral(double e1,double e2);
		
	/** Returns the sum all histogram values */
	double integral() const;
	/** Divide all energies by a constant factor */
	void divide(double div_,int option);
	void divide();
	/** Sets all histogram heights to zero */
	void reset();
		
	friend ostream& operator<<(ostream& s, Histogram & y);
	friend void CorrectSigma(Histogram &Sigma,double div_);        
	friend void MPI_VecReduce(Histogram &v1,Histogram &v2,int n,int type,int op,int rank,int comm);
      
private:
	
	double minE,maxE;
	int steps;
	double deltaE, *histE,*histDE;
          
};



#endif
