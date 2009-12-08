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

        jmatrix[0]=cos(dynVars.theta[site]);
	jmatrix[1]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
		-sin(dynVars.theta[site])*sin(dynVars.phi[site]));
	jmatrix[3]= -cos(dynVars.theta[site]);
	
	jmatrix[2]=conj(jmatrix[1]);
	
	
}


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)

{
	int	row, col, volume, i = 0, p;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j,k,dir;
	volume = ether.linSize;
	int level,gamma1,gamma2,spin1,spin2,index1, index2;
	vector<MatType> jmatrix(4,0);
	
	row=0;
	for (gamma1=0;gamma1<ether.numberOfOrbitals;gamma1++) {
		for (spin1=0;spin1<2;spin1++) {
			for (p = 0; p < volume; p++, row++) {
				index1 = spin1+gamma1*2;
				auxCreateJmatrix(jmatrix,dynVars,ether,p);
				
				matrix->rowptr[row] = i;
				matrix->colind[i] = row;		
				tmp =real(ether.J[gamma1]*jmatrix[spin1+2*spin1])*ether.modulus[p]+ether.potential[p+gamma1*volume];
				matrix->values[i] = (tmp-b)/a; /* -b because it's diagonal */
				i++;
				
				for (j = 0; j <  geometry.z(p); j++) {	
					k = geometry.neighbor(p,j);
					dir = geometry.scalarDirection(p,k,j);
					for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
						index2=spin1+2*gamma2; // spin2==spin1 because hoppings are diagonal in spin	
						matrix->colind[i]=k+index2*volume;
						matrix->values[i]=ether.bandHoppings[gamma1+gamma2*2+4*dir]/a;
						i++;
					}
				}
			
				gamma2=gamma1; // because J matrix is diagonal in orbital
				for (spin2=0;spin2<2;spin2++) {
					index2=spin2+2*gamma2;
					if (index1==index2) continue;
					matrix->colind[i]=p + index2*volume;
					matrix->values[i]=ether.J[gamma1]*jmatrix[spin1+2*spin2]*double(ether.modulus[p])/a;
					i++;
				}
			
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
	
	for (shift=0;shift<4;shift++) support.push_back(i+shift*volume);
	
	
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	int i;
	
	ether.typeofmodel="MODEL_KONDO_DMS_CUBIC";
	ether.hilbertSize=4*n;
	ether.nonzero= 4*(2*geometry.z(0)+2)*n;
	
	if (ether.numberOfOrbitals!=2) {
		cerr<<"I was expecting numberoforbitals=2 but instead you input ether.numberOfOrbitals="<<ether.numberOfOrbitals<<endl;
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	vectorPrint(ether.bandHoppings,"hoppings=",cerr);
	ether.energy1= -2*ether.D*maxElement(ether.bandHoppings)-maxElement(ether.J);
	ether.energy2= -ether.energy1;
	// calculate and add potential shift
	
	ether.energy1 += maxElement(ether.potential);
	ether.energy2 += maxElement(ether.potential);	
	aux.varTpem_a = 0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	
	
}
