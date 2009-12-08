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


#include "common.h"
	
void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
                 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)
{
        int     row, col, volume, i = 0, p;
        tpem_t hopping,hopping2;
        double  a=aux.varTpem_a, b=aux.varTpem_b, tmp,tmp2;
	int	j,iTmp;
	volume = ether.linSize;
	tpem_t S_ij;
	
	
	
	for (p = 0, row = 0; p < volume; p++, row++) {
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;		/* diagonal element */
		tmp = (ether.phononLambda*electronPhononTerm(p,geometry,dynVars,ether)+ether.potential[p]-b)/a;
		matrix->values[i] = tmp;
		
		i++;
		
		//cout<<"p="<<p<<endl;
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j);;
			if (iTmp>=0) hopping = ether.hoppings[iTmp]/a;
			else hopping= -1.0/a;

			col = geometry.neighbor(p,j);
			

			tmp=cos(0.5*dynVars.theta[p])*cos(0.5*dynVars.theta[col]);
			tmp2=sin(0.5*dynVars.theta[p])*sin(0.5*dynVars.theta[col]);
			 S_ij=tpem_t(tmp+tmp2*cos(dynVars.phi[p]-dynVars.phi[col]),
		 		-tmp2*sin(dynVars.phi[p]-dynVars.phi[col]));
			hopping2=hopping;
			if (p>col) {
				hopping2 = conj(hopping);
				
			} 
#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS_EX
			if (ether.phononDelta==0) {
				cerr<<"You set phononDelta=0.\n";
				cerr<<"This is not correct, ";
				cerr<<"please check cond-mat/0311200\n";
				cerr<<"FATAL: AT THIS POINT";
				cerr<<" "<<__FILE__<<" "<<__LINE__<<endl;
				exit(1);
			}
			tmp =  (1.0/ether.phononDelta+ether.phononAlpha*(
			electronPhononTerm(p,geometry,dynVars,ether)+
			electronPhononTerm(col,geometry,dynVars,ether)));
			hopping = tmp;
			hopping2 = hopping2 * hopping;
#endif	
		        matrix->colind[i] = col;
			matrix->values[i] = S_ij * hopping2;
			i++;
			
		}
	}
	if (i!=ether.nonzero) {
		cerr<<"i="<<i<<" but nonzero="<<ether.nonzero<<endl;
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}

	
}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{
	int j,col;
	support.push_back(i);
	for (j = 0; j < geometry.z(i); j++) {
		col = geometry.neighbor(i,j);
		support.push_back(col);
	}
	
}

void setHilbertParams(Parameters &ether, Aux &aux, Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	
	ether.typeofmodel="MODEL_KONDO_INF_ONEBAND_PHONONS";
	ether.hilbertSize=n;
	ether.nonzero=  (2*d + 1) * n;
	ether.energy1= -2*d-4*ether.phononLambda;
	ether.energy2= -ether.energy1;
	aux.varTpem_a = 0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	ether.classFieldList.push_back(1);
	
	
}	
