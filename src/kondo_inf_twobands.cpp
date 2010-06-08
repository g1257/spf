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

extern double calcPhonon(int ind,DynVars const &dynVars,Geometry const &g,Parameters const &ether,int what);


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
                 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)
  // Modified by IS Nov-08-04

{
	int	row, col, volume, i = 0, p,dir;
	tpem_t hopping,hopping2,bandHop;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp,tmp2;
	int	j,iTmp;
	volume = ether.linSize;
	tpem_t S_ij;
	vector<double> phonon_q1(volume);
	vector<double> phonon_q2(volume);
	vector<double> phonon_q3(volume);
	
	if (ether.isSet("adjusttpembounds")) {
		//cerr<<"setting a=1 and b=0\n";
		a=1.0; b=0.0;
	}
	
	Phonons<Parameters,Geometry> phonons(ether,geometry);
	
	for (p = 0, row = 0; p < volume; p++, row++) {
		
		
		phonon_q1[p]=phonons.calcPhonon(p,dynVars,0);
		phonon_q2[p]=phonons.calcPhonon(p,dynVars,1);
		phonon_q3[p]=phonons.calcPhonon(p,dynVars,2);	
		
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;		/* element aa */ 
		tmp = (ether.phononEjt[0]*phonon_q1[p]+ether.phononEjt[2]*phonon_q3[p]+ether.potential[p]-b)/a;
		matrix->values[i] = tmp;
		i++;
		
		matrix->colind[i] = row + volume;	/*  element ab */ 
		tmp = (ether.phononEjt[1]*phonon_q2[p])/a;
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
				
			if (p>col) hopping2 = conj(hopping);
			else hopping2=hopping;
			dir = int(j/2);
			
			bandHop=ether.bandHoppings[0+0*2+dir*4];
			hopping=hopping2 * bandHop ;		      
		        matrix->colind[i] = col;
			matrix->values[i] = hopping * S_ij;
			i++;
			
			bandHop=ether.bandHoppings[0+1*2+dir*4];
			hopping= hopping2 * bandHop;
		        matrix->colind[i] = col+volume;
			matrix->values[i] = hopping * S_ij;
			i++;
		}
	}
	for (p = 0; p < volume; p++, row++) {
	 	
		matrix->rowptr[row] = i;
		
		matrix->colind[i] = row ;	/* element bb */
		tmp = (-ether.phononEjt[2]*phonon_q3[p]+ether.phononEjt[0]*phonon_q1[p]+ether.potential[p]-b)/a;
		matrix->values[i] = tmp;
		i++;
		
		matrix->colind[i] = row - volume;		/*  element ba */
		tmp = (ether.phononEjt[1]*phonon_q2[p])/a;
		matrix->values[i] = tmp;
		i++;
		
		//cout<<"p="<<p<<endl;
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j);
			if (iTmp>=0) hopping = ether.hoppings[iTmp]/a;
			else hopping= -1.0/a;
			col = geometry.neighbor(p,j);
			tmp=cos(0.5*dynVars.theta[p])*cos(0.5*dynVars.theta[col]);
			tmp2=sin(0.5*dynVars.theta[p])*sin(0.5*dynVars.theta[col]);
			 S_ij=tpem_t(tmp+tmp2*cos(dynVars.phi[p]-dynVars.phi[col]),
		 		-tmp2*sin(dynVars.phi[p]-dynVars.phi[col]));
				
			if (p>col) hopping2 = conj(hopping);
			else hopping2=hopping;
			dir = int(j/2);
			
			bandHop=ether.bandHoppings[1+0*2+dir*4];
			hopping= hopping2 * bandHop;
			matrix->colind[i] = col;
			matrix->values[i] = hopping * S_ij;
			i++;
			
			bandHop=ether.bandHoppings[1+1*2+dir*4];
			hopping=hopping2 * bandHop;
			matrix->colind[i] = col+volume;
			matrix->values[i] = hopping * S_ij;
			i++;
		}
	}

}

void createHamiltonian(Geometry const &geometry, DynVars const &dynVars,
                 MyMatrix<std::complex<double> >& matrix,Parameters const &ether,Aux &aux,int type)
  // Modified by IS Nov-08-04

{
	int	col, volume, i = 0, p,dir;
	tpem_t hopping,hopping2,bandHop;
	double	 tmp,tmp2;
	int	j,iTmp;
	volume = ether.linSize;
	tpem_t S_ij;
	//static vector<double> phonon_q1(volume);
	//static vector<double> phonon_q2(volume);
	//static vector<double> phonon_q3(volume);
	
	
	Phonons<Parameters,Geometry> phonons(ether,geometry);
	for (p = 0; p < matrix.getRank(); p++)
		for (col = 0; col < matrix.getRank(); col++)
			matrix(p,col)=0;
	
	for (p = 0; p < volume; p++) {
		double phonon_q1=phonons.calcPhonon(p,dynVars,0);
		double phonon_q2=phonons.calcPhonon(p,dynVars,1);
		double phonon_q3=phonons.calcPhonon(p,dynVars,2);	
		matrix(p,p) = ether.phononEjt[0]*phonon_q1+ether.phononEjt[2]*phonon_q3+ether.potential[p];
		matrix(p+volume,p+volume) = -ether.phononEjt[2]*phonon_q3+ether.phononEjt[0]*phonon_q1+ether.potential[p];
		matrix(p,p+volume) = (ether.phononEjt[1]*phonon_q2);
		matrix(p+volume,p) = conj(matrix(p,p+volume));
		
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j);;
			if (iTmp>=0) hopping = ether.hoppings[iTmp];
			else hopping= -1.0;
			col = geometry.neighbor(p,j);
			tmp=cos(0.5*dynVars.theta[p])*cos(0.5*dynVars.theta[col]);
			tmp2=sin(0.5*dynVars.theta[p])*sin(0.5*dynVars.theta[col]);
			S_ij=tpem_t(tmp+tmp2*cos(dynVars.phi[p]-dynVars.phi[col]),
		 		-tmp2*sin(dynVars.phi[p]-dynVars.phi[col]));
				
			if (p>col) hopping2 = conj(hopping);
			else hopping2=hopping;
			dir = int(j/2);
			
			bandHop=ether.bandHoppings[0+0*2+dir*4];
			hopping=hopping2 * bandHop ;
			matrix(p,col) = hopping * S_ij;
			matrix(col,p) = conj(hopping * S_ij);
			
			bandHop=ether.bandHoppings[0+1*2+dir*4];
			hopping= hopping2 * bandHop;
		        matrix(p, col+volume)=hopping * S_ij;
			matrix(col+volume,p)=conj(matrix(p,col+volume));
			
			//
			bandHop=ether.bandHoppings[1+0*2+dir*4];
			hopping= hopping2 * bandHop;
			matrix(p+volume,col) =  hopping * S_ij;
			matrix(col,p+volume) =  conj(matrix(p+volume,col));
			
			bandHop=ether.bandHoppings[1+1*2+dir*4];
			hopping=hopping2 * bandHop;
			matrix(p+volume,col+volume) =  hopping * S_ij;
			matrix(col+volume,p+volume) = conj(matrix(p+volume,col+volume));
		}
		
		for (j = 0; j < geometry.z(p,2); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j,2);
			hopping= -1.0;
			if (iTmp>=0) hopping = ether.hoppings[iTmp];
			
			col = geometry.neighbor(p,j,2);
			tmp=cos(0.5*dynVars.theta[p])*cos(0.5*dynVars.theta[col]);
			tmp2=sin(0.5*dynVars.theta[p])*sin(0.5*dynVars.theta[col]);
			S_ij=tpem_t(tmp+tmp2*cos(dynVars.phi[p]-dynVars.phi[col]),
		 		-tmp2*sin(dynVars.phi[p]-dynVars.phi[col]));
				
			if (p>col) hopping2 = conj(hopping);
			else hopping2=hopping;
			dir = int(j/2);
			
			bandHop=ether.bandHoppings[0+0*2+dir*4]*ether.tprime;
			hopping=hopping2 * bandHop ;
			matrix(p,col) = hopping * S_ij;
			matrix(col,p) = conj(hopping * S_ij);
			
			bandHop=ether.bandHoppings[0+1*2+dir*4]*ether.tprime;
			hopping= hopping2 * bandHop;
		        matrix(p, col+volume)=hopping * S_ij;
			matrix(col+volume,p)=conj(matrix(p,col+volume));
			
			//
			bandHop=ether.bandHoppings[1+0*2+dir*4]*ether.tprime;
			hopping= hopping2 * bandHop;
			matrix(p+volume,col) =  hopping * S_ij;
			matrix(col,p+volume) =  conj(matrix(p+volume,col));
			
			bandHop=ether.bandHoppings[1+1*2+dir*4]*ether.tprime;
			hopping=hopping2 * bandHop;
			matrix(p+volume,col+volume) =  hopping * S_ij;
			matrix(col+volume,p+volume) = conj(matrix(p+volume,col+volume));
		}
	}
	//if (!matrix.isHermitian()) throw std::runtime_error("I'm barking\n");
	

}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{
	int j,col;
	support.push_back(i);
	support.push_back(i+geometry.volume());
	for (j = 0; j < geometry.z(i); j++) {
		col = geometry.neighbor(i,j);
		support.push_back(col);
		support.push_back(col+geometry.volume());
	}
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	
	ether.typeofmodel="MODEL_KONDO_INF_TWOBANDS";
	ether.hilbertSize=2*n;
	ether.nonzero= 4 * (2*d + 1) * n;
	if (!ether.isSet("spectrumbounds")) {
		ether.energy1= -2*d*maxElement(ether.bandHoppings)-maxElement(ether.potential);
		ether.energy2= -ether.energy1;
	}
	aux.varTpem_a = 1; //0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0; //0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	ether.classFieldList.push_back(1);
	
	
}	 
