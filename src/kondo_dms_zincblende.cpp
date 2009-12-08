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
	
void auxCreateJmatrix(vector<MatType> &jmatrix, DynVars const &dynVars,Parameters const &ether,int site)
{
	int i,j;
	
        jmatrix[0]=0.5*cos(dynVars.theta[site]);
        jmatrix[1]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                -0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));
        jmatrix[2]=0.0;
        jmatrix[3]=0.0;
        jmatrix[4]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                 0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));
        jmatrix[5]=cos(dynVars.theta[site])/6.0;
        jmatrix[6]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/3.0,
                 -sin(dynVars.theta[site])*sin(dynVars.phi[site])/3.0);
        jmatrix[7]=0.0;
        jmatrix[8]=0.0;
        jmatrix[9]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/3.0,
                 sin(dynVars.theta[site])*sin(dynVars.phi[site])/3.0);
        jmatrix[10]= -jmatrix[5];
        jmatrix[11]= MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                -0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));
        jmatrix[12]=0.0;
        jmatrix[13]=0.0;
        jmatrix[14]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                 0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));
        jmatrix[15]= -jmatrix[0];

	for (i=0;i<4;i++) {
		for (j=i+1;j<4;j++) {
			jmatrix[i+j*4]=conj(jmatrix[j+i*4]);
		}
	}
	
	for (i=0;i<16;i++) {
		jmatrix[i] *= ether.J[0];
		//cerr<<"i "<<jmatrix[i]<<endl;
	}

}


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)

{
	int	row, col, volume, i = 0, p;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j,k,dir;
	volume = ether.linSize;
	int level,gamma1,gamma2;
	vector<MatType> jmatrix(16,0);
	
	row=0;
	for (gamma1=0;gamma1<ether.numberOfOrbitals;gamma1++) {
		for (p = 0; p < volume; p++, row++) {
			
			auxCreateJmatrix(jmatrix,dynVars,ether,p);
			
			matrix->rowptr[row] = i;
			matrix->colind[i] = row;		
			tmp =real(jmatrix[gamma1+4*gamma1])*ether.modulus[p];
			matrix->values[i] = (tmp-b)/a; /* -b because it's diagonal */
			i++;
			
				
			for (j = 0; j <  geometry.z(p); j++) {	
				k = geometry.neighbor(p,j);
				dir = geometry.scalarDirection(p,k);
				for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
						
					matrix->colind[i]=k+gamma2*volume;
					matrix->values[i]=ether.bandHoppings[gamma1+gamma2*4+16*dir]/a;
					i++;
				}
			}
			if (p>=volume/2) continue;
			for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
				if (gamma1==gamma2) continue;
				matrix->colind[i]=p + gamma2*volume;
				matrix->values[i]=jmatrix[gamma1+4*gamma2]*double(ether.modulus[p])/a;
				i++;
			}
			
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
	int volume=geometry.volume();
	int shift;
	
	support.push_back(i);
	support.push_back(i+volume);
	support.push_back(i+2*volume);
	support.push_back(i+3*volume);
	
	
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	int i;
	
	ether.typeofmodel="MODEL_KONDO_DMS_ZINCBLENDE";
	ether.hilbertSize=4*n;
	ether.nonzero= ether.numberOfOrbitals*(1+ether.numberOfOrbitals*geometry.z(0)) * n+
		ether.numberOfOrbitals*n*0.5*(ether.numberOfOrbitals-1);
	if (!(geometry.latticeName()=="zincblende")) {
		cerr<<"I was expecting zincblende geometry but instead you input "<<geometry.latticeName()<<endl;
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	
	cerr<<ether.bandHoppings.size()<<endl;
	vectorPrint(ether.bandHoppings,"hoppings=",cerr);
	ether.energy1= -3-maxElement(ether.J);
	ether.energy2= -ether.energy1;
	aux.varTpem_a = 0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	
	
}
