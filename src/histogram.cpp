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

#include "histogram.h"
		 		 
using namespace std;
		 
	Histogram::~Histogram() 
	{

		// don't forget to deallocate all memory
		if (steps>0) {
			// cout<<"Histogram: Deallocating\n";
			delete [] histDE;
			delete [] histE;
		}
	}


	void Histogram::init(double minE_,double maxE_,int steps_)
	{
 		minE=minE_;
		maxE=maxE_;
		steps=steps_;
         
         
		// Some checking
    	 
		if (steps<=0 || minE>=maxE) {
			// handleError(__FILE__,__LINE__,E_HISTCON,"Cannot construct histogram");
			cerr<<"steps="<<steps<<" and minE="<<minE<<" and maxE="<<maxE<<endl;
			steps=0; // so we remember NOT to deallocate memory, see the destructor
			cerr<<"Histogram: "<<__FILE__<<":"<<__LINE__<<endl;
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

	int Histogram::add(double EValue) {
		int n;
		n=int(steps*(EValue-minE));
		n=int(n/(maxE-minE));
		// Dont remove this checking it's very important!
		if (n>=steps || n<0) return 1;	  
		histDE[n]++;
	  	return 0;
	}

			  
	int Histogram::add(double xValue,double EValue) {
		int n;
		n=int(steps*(xValue-minE));
		n=int(n/(maxE-minE));	
		// Dont remove this checking it's very important!
		if (n>=steps || n<0) return 1;	    
		histDE[n]+=EValue;
	  	return 0;	  		  
	} 


	// sum all histogram values from energy e1 to energy e2
	double Histogram::integral(double e1,double e2)  const
	{
		int i;
		double s=0;
			
		for (i=0;i<steps;i++) {
			if (histE[i]>=e1)  s+=histDE[i];
			if (histE[i]>=e2) break;
		}
		return s;
	}	  	  
	  
	// sum all histogram values
	double Histogram::integral()  const
	{
		int i;
		double s=0;
			
		for (i=1;i<steps;i++) s+=histDE[i];
		return s;
	}
		
	// divide all energies by a constant factor
	void Histogram::divide(double div_,int option) 
	{
		int i;
		if (option) divide();
			
		if (div_<=0) {
			cerr<<"\nERROR AT: ";
			cerr<<__FILE__<<" : "<<__LINE__<<" Internal, div_="<<div_<<endl;
			exit(1);
		}
			
		for (i=0;i<steps;i++) histDE[i] = (double)histDE[i]/div_;

	}
	void Histogram::divide() 
	{
		int i;
			
		double div_=(double)(maxE-minE)/steps;
			
		for (i=0;i<steps;i++) histDE[i] = (double)histDE[i]/div_;

	}	  
		
	void Histogram::reset()
	{
			
		for (int n=0;n<steps+1;n++) {
			histE[n]=minE+n*deltaE;
			histDE[n]=0;
		}
	}	    
	
      
