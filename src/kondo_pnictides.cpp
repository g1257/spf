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


void createHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 MyMatrix<std::complex<double> >& matrix,Parameters const &ether,Aux &aux,int type)
{
	size_t volume = ether.linSize;
	size_t norb = ether.numberOfOrbitals;
	size_t dof = norb * 2; // the 2 comes because of the spin
	vector<MatType> jmatrix(dof*dof);
	
	for (size_t gamma1=0;gamma1<matrix.getRank();gamma1++) 
		for (size_t p = 0; p < matrix.getRank(); p++) 
			matrix(gamma1,p)=0;
	
	for (size_t gamma1=0;gamma1<dof;gamma1++) {
		for (size_t p = 0; p < volume; p++) {
			
			auxCreateJmatrix(jmatrix,dynVars,ether,p);

			size_t spin1 = size_t(gamma1/2);
			size_t orb1 = gamma1 % 2;
			matrix(p+gamma1*volume,p+gamma1*volume) = real(jmatrix[spin1+2*spin1]) + ether.potential[p];
				
			for (int j = 0; j <  geometry.z(p); j++) {	
				if (j%2!=0) continue;	
				int k = geometry.neighbor(p,j);
				
				int dir = int(j/2);
				for (size_t orb2=0;orb2<norb;orb2++) {
					size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
					matrix(p+gamma1*volume,k+gamma2*volume)=
						ether.bandHoppings[orb1+orb2*norb+norb*norb*dir];
					matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
				}
                     	}
			//if (geometry.z(p,2)!=4 || geometry.z(p)!=4) throw std::runtime_error("neighbours wrong\n");
			for (int j = 0; j <  geometry.z(p,2); j++) {
				if (j%2!=0) continue;	
				int k = geometry.neighbor(p,j,2);
				int dir = int(j/2)+2;
				//std::cerr<<"Neigbors "<<p<<" "<<k<<"\n";
				for (size_t orb2=0;orb2<norb;orb2++) {
					size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
					matrix(p+gamma1*volume,k+gamma2*volume)=
						ether.bandHoppings[orb1+orb2*norb+norb*norb*dir];
					matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
				}
                     	}
			
			for (size_t spin2=0;spin2<2;spin2++) {
				if (spin1==spin2) continue; // diagonal term already taken into account
				size_t gamma2 = orb1+spin2*norb; // orb2 == orb1 here
				matrix(p+gamma1*volume,p + gamma2*volume)=jmatrix[spin1+2*spin2];
				matrix(p + gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,p + gamma2*volume));
			}
		}
	}
// 	std::cerr<<matrix;
// 	throw std::runtime_error("testing\n");
	
	if (!matrix.isHermitian()) {
		
		throw std::runtime_error("Problem matrix non hermitian!!\n");
	}
}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{
	int volume=geometry.volume();
	
	support.push_back(i);
	support.push_back(i+volume);
	support.push_back(i+2*volume);
	support.push_back(i+3*volume);
        
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume();
	
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
	aux.varTpem_a = 1; //0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0; //0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	
	
}
