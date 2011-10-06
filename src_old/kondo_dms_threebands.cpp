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

	MatType a1,a2;
	
	a1=cos(dynVars.theta[site]);
	a2=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
		-sin(dynVars.theta[site])*sin(dynVars.phi[site]));
	for (i=0;i<3;i++) jmatrix[7*i]=a1;
	for (i=3;i<6;i++) jmatrix[7*i]= -a1;
	for (i=0;i<3;i++) jmatrix[7*i+3]=a2;	
	for (i=3;i<6;i++) jmatrix[7*i-3]= conj(a2);
	for (i=0;i<36;i++) {
		jmatrix[i] *= ether.J[0];
		//cerr<<"i "<<jmatrix[i]<<endl;
	}

}


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)

{
	int	row, col, volume, i = 0, p;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j,k,dir,dir2;
	volume = ether.linSize;
	int level,gamma1,gamma2;
	MatType temp;
	vector<MatType> jmatrix(36,0);
	
	row=0;
	for (gamma1=0;gamma1<ether.numberOfOrbitals;gamma1++) {
		for (p = 0; p < volume; p++, row++) {
			
			auxCreateJmatrix(jmatrix,dynVars,ether,p);
			
			matrix->rowptr[row] = i;
			matrix->colind[i] = row;		
			tmp =real(jmatrix[gamma1+ether.numberOfOrbitals*gamma1])*double(ether.modulus[p]);
			tmp =0;
			matrix->values[i] = (tmp-b)/a; /* -b because it's diagonal */
			i++;
			
				
			for (j = 0; j <  geometry.z(p); j++) {	
				k = geometry.neighbor(p,j);
				dir = geometry.scalarDirection(p,k);
				for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
						
					matrix->colind[i]=k+gamma2*volume;
					dir2 = dir;
					if (dir2 >3) dir2 -= 4;
					
					temp = ether.bandHoppings[gamma1+gamma2*ether.numberOfOrbitals+
	ether.numberOfOrbitals*ether.numberOfOrbitals*dir2]/a;
					if (dir>3) temp = conj(temp);
					
					matrix->values[i]=temp;
					i++;
				}
			}
			if (p>=volume/2) continue;
			for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
				if (gamma1==gamma2) continue;
				matrix->colind[i]=p + gamma2*volume;
				matrix->values[i]=jmatrix[gamma1+ether.numberOfOrbitals*gamma2]*double(ether.modulus[p])/a;
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
	
	ether.typeofmodel="MODEL_KONDO_DMS_SIXBANDS";
	ether.hilbertSize=6*n;
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
	ether.nOfFields=1;
	
	
}
