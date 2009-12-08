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
        double  a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j;
	
	volume = ether.linSize;
	if (ether.isSet("adjusttpembounds")) {
		a=1.0; b=0.0;
	}
	
	for (p = 0, row = 0; p < volume; p++, row++) {
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;
		tmp = (ether.potential[p]- b) / a;		/* diagonal element */				
		matrix->values[i] = tmp;
		i++;
		
		matrix->colind[i] = row + volume;	/* off-diag element (delta-term onsite) */
		tmp=0;
		if (ether.isSet("bcsonsite")) tmp = ether.bcsV[p] * dynVars.bcsDelta[p]/a;
		matrix->values[i] = MatType(tmp*cos(dynVars.bcsPhi[p][0]),tmp*sin(dynVars.bcsPhi[p][0]));
		i++;
		if (!ether.isSet("bcsonsite")) {
			for (j = 0; j < geometry.z(p); j++) {	/* off-diag element (delta-term nn) */
				
				col = geometry.neighbor(p,j);
				if (j==0 || j %2 ==0) {
					tmp = dynVars.bcsPhi[p][j/2]; // forward link
					hopping = ether.bcsV[p+j*volume/2] * dynVars.bcsDelta[p]/a;
				} else {
					tmp = dynVars.bcsPhi[col][(j-1)/2]; //backward link
					hopping = ether.bcsV[col+(j-1)*volume/2] * dynVars.bcsDelta[p]/a;
				}
		        	matrix->colind[i] = col+volume;
				hopping2 = tpem_t(cos(tmp),sin(tmp));
				matrix->values[i] = hopping*hopping2;			
				i++;
			}
		}
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			hopping = -1.0/a;
			col = geometry.neighbor(p,j);
			if (p>col) hopping = conj(hopping);
		        matrix->colind[i] = col;
			matrix->values[i] = hopping;			
			i++;
		}
	}
	for (p = 0; p < volume; p++, row++) {
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;		/* diagonal element */
		tmp = (ether.potential[p]- b) / a;
		matrix->values[i] = tmp;
		i++;
		matrix->colind[i] = row - volume;	/* off-diag element (delta term) */
		tmp=0;
		if (ether.isSet("bcsonsite")) tmp = ether.bcsV[p] * dynVars.bcsDelta[p]/a;
		matrix->values[i] = tpem_t(tmp*cos(dynVars.bcsPhi[p][0]), -tmp*sin(dynVars.bcsPhi[p][0]));
		i++;
		if(!ether.isSet("bcsonsite")) {
			for (j = 0; j < geometry.z(p); j++) {	/* off-diag element (delta-term nn) */
				
				col = geometry.neighbor(p,j);
				if (j==0 || j %2 ==0) {
					tmp = dynVars.bcsPhi[p][j/2]; // forward link
					hopping = ether.bcsV[p+j*volume/2] * dynVars.bcsDelta[p]/a;
				} else {
					tmp = dynVars.bcsPhi[col][(j-1)/2]; //backward link
					hopping = ether.bcsV[p+(j-1)*volume/2] * dynVars.bcsDelta[p]/a;
				}
		        	matrix->colind[i] = col;
				hopping2 = tpem_t(cos(tmp),-sin(tmp));
				matrix->values[i] = hopping*hopping2;			
				i++;
			}
		}	
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			hopping = 1.0/a;
			
			col = geometry.neighbor(p,j) + volume;
			if (p>col-volume) hopping = conj(hopping);
			matrix->colind[i] = col;
			matrix->values[i] = hopping;
			i++;
		}
	}
	if (i!=ether.nonzero) {
		cerr<<"i="<<i<<" but nonzero="<<ether.nonzero<<endl;
		exit(1);
	}
	
}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{
	int j;
	
	support.push_back(i);
	support.push_back(i + geometry.volume());
	
	for (j = 0; j < geometry.z(i); j++) {
		support.push_back(geometry.neighbor(i,j));
		support.push_back(geometry.neighbor(i,j)+geometry.volume());
	}
	
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume();
	
	ether.typeofmodel="MODEL_KONDO_BCS";
	ether.hilbertSize=2*n;
	if (ether.isSet("bcsonsite")) {
		ether.nonzero= 2 * (geometry.z(0) + 2) * n;
	} else {
		ether.nonzero= 2 * (2*geometry.z(0) + 2) * n;
	}
	if (!ether.isSet("spectrumbounds")) {
		ether.energy1 = -sqrt(square(2*geometry.dim()+fabs(ether.potential[0]))+square(ether.bcsDelta0));
		ether.energy2 = -ether.energy1;
	}
	aux.varTpem_a = 0.5*(ether.energy2 - ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2 + ether.energy1);
	ether.classFieldList.push_back(2);
	
	
}


