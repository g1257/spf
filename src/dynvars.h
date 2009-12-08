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




#ifndef DYNVARS_H_
#define DYNVARS_H_

/** \brief Dynamical Variables for the Monte Carlo simulation 

* This is a public class with only member data. Don't add member functions.
* Don't add pointers to avoid having to write copy constructors or
* assignment operators. */
struct DynVars {
	
	
	std::vector<double> theta,phi;
	std::vector<std::vector<double> > phonons;
	std::vector<double> bcsDelta;
	std::vector<vector<double> > bcsPhi;
	
	std::vector<double> Tx,Ty,Tz; //orbital pseudospin
};

inline void copyDynVars(DynVars &dv1,int i,DynVars const &dv2,int j) 
{
	if (dv1.theta.size()>i && dv2.theta.size()>j) {
		dv1.theta[i]=dv2.theta[j];
		dv1.phi[i]=dv2.phi[j];
	}
	if (dv1.phonons.size()>i && dv2.phonons.size()>j) 
		dv1.phonons[i]=dv2.phonons[j];
	if (dv1.bcsDelta.size()>i && dv2.bcsDelta.size()>j) 
		dv1.bcsDelta[i]=dv2.bcsDelta[j];
	if (dv1.bcsPhi.size()>i && dv2.bcsPhi.size()>j) 
		dv1.bcsPhi[i]=dv2.bcsPhi[j];
}

inline void copyDynVars(DynVars &dv1,DynVars const &dv2)
{
	if (dv2.theta.size()>0) dv1.theta=dv2.theta;
	if (dv2.phi.size()>0)  dv1.phi=dv2.phi;
	if (dv2.phonons.size()>0) dv1.phonons=dv2.phonons;
	if (dv2.bcsDelta.size()>0) dv1.bcsDelta=dv2.bcsDelta;
	if (dv2.bcsPhi.size()>0) dv1.bcsPhi=dv2.bcsPhi;
} 
#endif

