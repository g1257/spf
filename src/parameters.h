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




#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <iostream>
#include <vector>
#include "basic_types.h"

using std::vector;
using std::string;
using std::ostream;

const std::string ansiRed="\033[0;31;40m";
const std::string ansiGreen="\033[0;32;40m";
const std::string ansiBlue="\033[0;34;40m";
const std::string ansiMagenta="\033[0;35;40m";
const std::string ansiCyan="\033[0;36;40m";
const std::string ansiWhite="\033[0;37;40m";
const std::string ansiYellow="\033[1;33;40m";
const std::string ansiReset="\033[0m";
/* [1; for bold [3; for italics */

/** \brief A class to store parameters 

* If parameters change frequently store them in the class Aux.
* All member data is public. Don't add pointer data members to
* avoid having to write copy constructor and assignment operator */

struct Parameters {
	int iterTherm,iterEffective,iterUnmeasured;
	double window,beta;
	size_t numberOfBetas;
	std::vector<double> betaVector;
	int D,linSize,mcflag,startType,startLevel,conc;
	vector<double> JafVector;
	double jafCenter,jafDelta;
	size_t numberOfJafConfigs;
	vector<double> potential;
	std::string potentialFile; //< the filename from where to load the potential
	std::string jafFile; //< the filename from where to load the jaf vector
	std::string rootname;
	std::string version,typeofmodel;
	int histSteps,carriers,hilbertSize,nonzero;
	/** This classFieldList is a list of indices that denote different classical Fields
	 ** 0 is O(3) spin, 1 is phonons, 2 is bcs delta/phi */
	std::vector<int> classFieldList;
	double energy1, energy2, phononAlpha, phononDelta, phononGj,maxPhonons,hist1,hist2;
	int numberOfOrbitals;
	 
	vector<MatType> hoppings;
	vector<MatType> bandHoppings;
	vector<int> modulus;
	std::string options;
	std::string dynVarsInputFile;
	vector<double> J,phononEjt,phononEd;
	double tprime; //next nearest neighbour hopping
	double jafprime; // next nearest neighbor direct exchange coupling
	double tsecond; // next next nearest neighbor hopping
	vector<double> magnetic; // magnetic field
	vector<double> bcsV;
	double bcsDelta0;
	std::string typeoflattice; //> type of lattice
	
	double tpem_epsProd,tpem_epsTrace;
	int tpem_cutoff,tpem;
	int mpiRank,mpiTpemRank,mpiNop1,mpiNop2,mpiSize;
	std::vector<int> localRank;
	std::vector<size_t> localSize;
	
	void print(ostream &s,char prefix='#');
	bool isSet(char const *what) const;
	void setOption(char const *what);
	void setVersion(char s[]);
	
	void welcome();
	Parameters();
#ifdef MODEL_KONDO_PHONONS
	double phononLambda;
#endif				
		
	
};

#endif

