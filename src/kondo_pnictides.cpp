/*******************************************************************
 **         THIS FILE IS PART OF THE SPF (SPIN-PHONON-FERMION)    ** 
 **                          COMPUTER PROGRAM                     **
 **                                                               **
 ********************** SPF VERSION 6.2. ***************************
 **                                                               **
 **       For more info see: http://mri-fre.ornl.gov/spf          **
 ******************************************************************/

#include "common.h"
	
void auxCreateJmatrix(vector<MatType> &jmatrix,DynVars const &dynVars,Parameters const &ether,int site)
{
	
        jmatrix[0]=0.5*cos(dynVars.theta[site]);

        jmatrix[1]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
                sin(dynVars.theta[site])*sin(dynVars.phi[site]));     
   
        jmatrix[2]=conj(jmatrix[1]);

        jmatrix[3]= -cos(dynVars.theta[site]);
 
        for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= ether.J[0]; 
}


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)
{
	int	row, col, volume, i = 0, p;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j,k,dir;
	volume = ether.linSize;
	int level,gamma1,gamma2;
	int dof = ether.numberOfOrbitals * 2; // the 2 comes because of the spin
	vector<MatType> jmatrix(dof*dof);
        
	  if (ether.isSet("adjusttpembounds")) {
                a=1.0; b=0.0;
        }

	row=0;
	
	for (gamma1=0;gamma1<dof;gamma1++) {
			for (p = 0; p < volume; p++, row++) {
			
			auxCreateJmatrix(jmatrix,dynVars,ether,p);

			size_t spin1 = size_t(gamma1/2);
			size_t orb1 = gamma1 % 2;
			matrix->rowptr[row] = i;
			matrix->colind[i] = row;		
			tmp =real(jmatrix[spin1+2*spin1]) + ether.potential[p];
//                        tmp =real(jmatrix[gamma1+6*gamma1])*double(ether.modulus[p]) + ether.potential[p];
  			
			matrix->values[i] = (tmp-b)/a; /* -b because it's diagonal */
			i++;
			
				
			for (j = 0; j <  geometry.z(p); j++) {	
				k = geometry.neighbor(p,j);
				dir = geometry.scalarDirection(p,k,j);
//                               for (gamma1=0;gamma1<ether.numberOfOrbitals;gamma1++) {
				for (gamma2=0;gamma2<dof;gamma2++) {
					matrix->colind[i]=k+gamma2*volume;
					matrix->values[i]=ether.bandHoppings[gamma1+gamma2*dof+dof*dof*dir]/a;
					i++;
				}
//			}
			
                     }
			
			for (gamma2=0;gamma2<dof;gamma2++) {
				size_t spin2 = size_t(gamma2/2);
				size_t orb2 = gamma2 % 2;
				if (orb1!=orb2) continue; // diagonal in orbitals
				if (spin1==spin2) continue; // diagonal term already taken into account 
				matrix->colind[i]=p + gamma2*volume;
				matrix->values[i]=jmatrix[spin1+2*spin2]/a;
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
	
	ether.typeofmodel="MODEL_KONDO_PNICTIDES";
	ether.hilbertSize=2*n*ether.numberOfOrbitals;
	ether.nonzero= square(ether.numberOfOrbitals)*(1+geometry.z(0)) * n;
	
	
//	cerr<<ether.bandHoppings.size()<<endl;
//	vectorPrint(ether.bandHoppings,"hoppings=",cerr);
//	ether.energy1= -10-maxElement(ether.J);
	if (!ether.isSet("spectrumbounds")) {
        	ether.energy1= -18.5;
		ether.energy2= 20.5;
	}
	aux.varTpem_a = 0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	
	
}
