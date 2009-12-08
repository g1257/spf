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
	int	j,iTmp;

	volume = ether.linSize;
	if (ether.isSet("adjusttpembounds")) {
		a=1.0; b=0.0;
	}
	
	for (p = 0, row = 0; p < volume; p++, row++) {
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;
		tmp=0;		/* diagonal element */
#ifdef MODEL_KONDO_FINITE_PHONONS
		tmp= ether.phononLambda*electronPhononTerm(p,geometry,dynVars,ether);
#endif				
		tmp = (tmp-ether.J[type]  * cos (dynVars.theta[p]) * ether.modulus[p] + ether.potential[p]- b) / a;
		matrix->values[i] = tmp;
		i++;
		matrix->colind[i] = row + volume;	/* off-diag element */
		tmp = -ether.J[type] * sin (dynVars.theta[p]) * ether.modulus[p];
		tmp /= a;
		matrix->values[i] = tpem_t(tmp * cos (dynVars.phi[p]),
						tmp * sin (dynVars.phi[p]));
		i++;
		//cout<<"p="<<p<<endl;
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j);
			if (iTmp>=0) hopping2 = ether.hoppings[iTmp]/a;
			else hopping2 = -1.0/a;
			hopping = hopping2;
			col = geometry.neighbor(p,j);
			if (p>col) hopping = conj(hopping);
#ifdef MODEL_KONDO_FINITE_PHONONS			
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
			hopping2 = tmp;
			hopping = hopping2 * hopping;
#endif			
		        matrix->colind[i] = col;
			matrix->values[i] = hopping;
			//cout<<"  -->"<<col<<" ";
			//cout<<"("<<tpem_real(hopping)<<","<<tpem_imag(hopping)<<")\n";
			
			i++;
			if (ether.isSet("tprime")) {
				/* next nearest neighbor hopping */
				col=geometry.neighbor(p,j);
				hopping = -ether.tprime/a;
				if (p>col-volume) hopping = conj(hopping);
				matrix->colind[i] = col;
				matrix->values[i] = hopping;
				//cout<<"("<<tpem_real(hopping)<<","<<tpem_imag(hopping)<<") i="<<i<<" nonzero="<<ether.nonzero<<"\n";
				i++;
			} // tprime
				if (ether.isSet("tsecond")) {
				/* next nearest neighbor hopping */
				col=geometry.neighbor(p,j,2);
				hopping = -ether.tsecond/a;
				if (p>col-volume) hopping = conj(hopping);
				matrix->colind[i] = col;
				matrix->values[i] = hopping;
				//cout<<"("<<tpem_real(hopping)<<","<<tpem_imag(hopping)<<") i="<<i<<" nonzero="<<ether.nonzero<<"\n";
				i++;
			} // tsecond
		}
	}
	for (p = 0; p < volume; p++, row++) {
		matrix->rowptr[row] = i;
		matrix->colind[i] = row;		/* diagonal element */
		tmp=0;
#ifdef MODEL_KONDO_FINITE_PHONONS
		tmp= ether.phononLambda*electronPhononTerm(p,geometry,dynVars,ether);
#endif	
		tmp = (tmp+ether.J[type] * cos (dynVars.theta[p])* ether.modulus[p] + ether.potential[p]- b) / a;
		matrix->values[i] = tmp;
		i++;
		matrix->colind[i] = row - volume;	/* off-diag element */
		tmp = -ether.J[type] * sin (dynVars.theta[p]) * ether.modulus[p];
		tmp /= a;
		matrix->values[i] = tpem_t(tmp * cos (dynVars.phi[p]),
						-tmp * sin (dynVars.phi[p]));
		i++;
		//cout<<"p="<<p<<endl;
		for (j = 0; j < geometry.z(p); j++) {	/* hopping elements */
			iTmp=geometry.borderId(p,j);
			if (iTmp>=0) hopping2 = ether.hoppings[iTmp]/a;
			else hopping2 = -1.0/a;
			
			hopping = hopping2;
			//hopping= tpem_cast((-1.0-b)/a,0.0);
			col = geometry.neighbor(p,j) + volume;
			if (p>col-volume) hopping = conj(hopping);
			
#ifdef MODEL_KONDO_FINITE_PHONONS			
			tmp =  (1.0/ether.phononDelta+ether.phononAlpha*(
			electronPhononTerm(p,geometry,dynVars,ether)+
			electronPhononTerm(col-volume,geometry,dynVars,ether)));
			hopping2 = tmp;
			hopping = hopping2 * hopping;
#endif	

			matrix->colind[i] = col;
			matrix->values[i] = hopping;
			//cout<<"  -->"<<col<<" ";
			//tpem_print(hopping,cout);
			//cout<<endl;
			i++;
			if (ether.isSet("tprime")) {
				/* next nearest neighbor hopping */
				col=geometry.neighbor(p,j)+volume;
				hopping = -ether.tprime/a;
				if (p>col-volume) hopping = conj(hopping);
				matrix->colind[i] = col;
				matrix->values[i] = hopping;
				i++;
				//cout<<" down ("<<tpem_real(hopping)<<","<<tpem_imag(hopping)<<") i="<<i<<" nonzero="<<ether.nonzero<<"\n";
			} // tprime
			if (ether.isSet("tsecond")) {
				/* next next nearest neighbor hopping */
				col=geometry.neighbor(p,j)+volume;
				hopping = -ether.tsecond/a;
				if (p>col-volume) hopping = conj(hopping);
				matrix->colind[i] = col;
				matrix->values[i] = hopping;
				i++;
				//cout<<" down ("<<tpem_real(hopping)<<","<<tpem_imag(hopping)<<") i="<<i<<" nonzero="<<ether.nonzero<<"\n";
			} // tsecond
		}
	}
	if (i!=ether.nonzero) {
		cerr<<"i="<<i<<" but nonzero="<<ether.nonzero<<endl;
		exit(1);
	}
	
}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{

	support.push_back(i);
	support.push_back(i + geometry.volume());
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	
	ether.typeofmodel="MODEL_KONDO_FINITE";
#ifdef MODEL_KONDO_FINITE_PHONONS
	ether.verson="MODEL_KONDO_FINITE_PHONONS";
#endif
	ether.hilbertSize=2*n;
	ether.nonzero= 2 * (geometry.z(0) + 2) * n;
	
	if (ether.isSet("tsecond")) ether.nonzero += 4*d*n;
	if (ether.isSet("tprime")) {
		ether.nonzero += 4*d*n;	
		if (fabs(ether.tprime)>0.5 || ether.tprime>0) {
			if (ether.tpem>0) {
				cerr<<"TPEM AND TPRIME ONLY WITH -0.5<=TPRIME<=0\n";
				cerr<<"AT THIS POINT IN FILE : "<<__LINE__<<endl;
				exit(1);
			}
		}
	}
	ether.energy1= -2*ether.D-maxElement(ether.J)-maxElement(ether.potential)+2*d*fabs(ether.tprime);	
	ether.energy2= 2*ether.D+maxElement(ether.J)+maxElement(ether.potential) +2*d*fabs(ether.tprime);
	aux.varTpem_a = 0.5*(ether.energy2 - ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2 + ether.energy1);
	ether.classFieldList.push_back(0);
	if (ether.isSet("verbose") && ether.mpiRank==0) {
		cerr<<" potential="<<maxElement(ether.potential);
		cerr<<" tprime="<<ether.tprime<<endl;
		cerr<<" tsecond="<<ether.tsecond<<endl;
	}
#ifdef MODEL_KONDO_FINITE_PHONONS
	ether.classFieldList.push_back(1);
	ether.energy1 -= 4*ether.phononLambda;
	ether.energy2 += 4*ether.phononLambda;
	aux.varTpem_a += 4*ether.phononLambda;
#endif
	
	
}
