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




#ifndef COMMON_H_
#define COMMON_H_

#include "tpemplus.h"
#include "geometry.h"
#include "io.h"
#include <functional>
 int loadFromFile(char const *filename,Geometry &geometry,DynVars
&dynVars,Parameters &ether,Aux &aux);
 void setupVariables(Geometry const &geometry,DynVars &dynVars,Parameters
&ether,Aux &aux);
 void doMonteCarlo(Geometry const &geometry,DynVars &dynVars,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions);
 void doMeasurements(int iter,DynVars const &dynVars,Geometry const &geometry,Io &io,
		Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions);
 void diag(vector<double> &eigOneBand,Geometry const &geometry,
        DynVars const &dynVars,Parameters const &ether, Aux &aux,char jobz='N');
 void tmpValues(double &a,double &b,double &mu,double &beta,int option);
 double effective_action (double x);
void setupHamiltonian(MyMatrix<MatType> &matrix,Geometry const &geometry,DynVars const &dynVars,
		Parameters const &ether,Aux &aux,int type);
void calcMoments(DynVars const &dynVars,Geometry const &geometry,Parameters const
&ether,Aux &aux,TpemOptions const &tpemOptions);
void tpemOptionsFill(TpemOptions &tpemOptions,Parameters const &ether);
void customConfig(Parameters &ether,Aux &aux,TpemOptions const &tpemOptions,int nop3,double shift=0);
void tpemOptionsFill(TpemOptions &tpemOptions,Parameters const &ether);
int spf_entry(int argc,char *argv[],int mpiRank, int mpiSize);

#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
double electronPhononTerm(int p,Geometry const &geometry, DynVars const
&dynVars,Parameters const &ether);
#endif

// from adjustments.cpp
void tpem_sparse_scale(tpem_sparse *matrix,double a,double b);
int tpemAdjustBounds(tpem_sparse *matrix,Parameters const &ether, Aux &aux);
double nOfElectrons(double mu,double beta,vector<double> const &eig);
double nOfElectPrime(double mu,double beta,vector<double> const &eig);
void adjChemPot(vector<double> const &eig,Parameters const &ether,Aux &aux);
void adjChemPotTpem(Parameters const &ether, Aux &aux,TpemOptions const &tpemOptions);

// from groundstate.cpp
void calcGroundState(Geometry const &geometry,DynVars &dynVars,
Parameters const &ether,Aux &aux);

// from manybody.cpp
void accNOfOmega(Geometry const &geometry,DynVars const &dynVars, Parameters const
&ether,Aux &aux);
void accLdos(Geometry const &geometry,DynVars const &dynVars, Parameters const
&ether,Aux &aux);
void accAkw(Geometry const &geometry,DynVars const &dynVars, Parameters const
&ether,Aux &aux);
void accLcd(Geometry const &geometry,DynVars const &dynVars,Parameters const
&ether,Aux &aux);
double calcKinetic(DynVars const &dynVars,Geometry const &geometry,Parameters
const &ether,Aux &aux);
double measure_kinetic(Geometry const &geometry,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions);
void kTpemMomentsOptical(Geometry const &geometry, Parameters const &ether,Aux
&aux,TpemOptions const &tpemOptions);
void accOptical(Geometry const &geometry,DynVars const &dynVars,Parameters const
&ether,Aux &aux);
void kTpemMomentsCl(Geometry const &geometry, Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions);
double calcNumber(DynVars const &dynVars,Geometry const &geometry,Parameters const &ether,Aux &aux,int band);
void accChargeCorrelation(size_t g1,size_t g2,Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux);
void accOrbitalAngles(Geometry const &geometry,Parameters const &ether,Aux &aux);
void accOrbitalCorrelation(size_t g1,size_t g2,Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux);


#endif

