
/** \ingroup SPF */
/*@{*/

/*! \file ThreeOrbitalTerms.h
 *
 * Additional terms needed by the 3-orbital
 * Hamiltonian, such as t7, t8 and deltaxy
 *
 * author: Q.L.
 * need to add AUTHORS file at the top
 * need to ask all AUTHORS how they want to appear there
 * FIXME
 */

#ifndef THREE_ORBITAL_TERMS_H
#define THREE_ORBITAL_TERMS_H
#include "Vector.h"

namespace Spf {
	template<int norb_,typename MatrixType,typename ParametersModelType,
	typename GeometryType>
	class ThreeOrbitalTerms {
		typedef typename ParametersModelType::RealType RealType;
	public:
		ThreeOrbitalTerms(
			const ParametersModelType& mp,
			const GeometryType& geometry) {}
			void operator()(MatrixType& ham) const { }
			
			RealType hopping(size_t isite,size_t dir2,size_t iorb,size_t iorb2) const
			{
				return 0;
			}
	}; // class ThreeOrbitalTerms

	template<typename MatrixType,typename ParametersModelType,
		typename GeometryType>
	class ThreeOrbitalTerms<3,MatrixType,ParametersModelType,GeometryType> {
		static const size_t SPINS = 2; // 2 spins
		typedef typename MatrixType::value_type FieldType;
		typedef typename ParametersModelType::RealType RealType;
		enum {
			DIR_X=GeometryType::DIRX,
			DIR_Y=GeometryType::DIRY,
			DIR_XPY=GeometryType::DIRXPY,
			DIR_XMY=GeometryType::DIRXMY};

	public:
		ThreeOrbitalTerms(
				const ParametersModelType& mp,
				const GeometryType& geometry)
		: mp_(mp),geometry_(geometry)
		{}

		void operator()(MatrixType& ham) const
		{
			/*--------- hoppings t7 ~ t8 ------------*/
			size_t volume = geometry_.volume();
			for (size_t ispin=0;ispin<SPINS;ispin++)
				for (size_t isite = 0; isite <volume; isite++) {

					// xz -> xy, +x
					size_t iorb=0;
					size_t iorb2=2;
					size_t ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					size_t isite2=geometry_.getNeighbour(isite,DIR_X).first;
					size_t iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					FieldType hopping = signHop(isite)*mp_.t7;
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> xz, +x
					iorb=2; iorb2=0;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_X).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*mp_.t7;
					updateHamiltonian(ham,ix,iy,hopping);

					// yz -> xy, +y
					iorb=1; iorb2=2;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_Y).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*mp_.t7;
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> yz, +y
					iorb=2; iorb2=1;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_Y).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*mp_.t7;
					updateHamiltonian(ham,ix,iy,hopping);

					// xz -> xy, +x+y
					iorb=0; iorb2=2;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XPY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*mp_.t8;
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> xz, +x+y
					iorb=2; iorb2=0;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XPY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(-mp_.t8);
					updateHamiltonian(ham,ix,iy,hopping);

					// xz -> xy, +x-y
					iorb=0; iorb2=2;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XMY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(-mp_.t8); // ***********************
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> xz, +x-y
					iorb=2; iorb2=0;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XMY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(mp_.t8); // ***********************
					updateHamiltonian(ham,ix,iy,hopping);

					// yz -> xy, +x+y
					iorb=1; iorb2=2;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XPY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*mp_.t8;
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> yz, +x+y
					iorb=2; iorb2=1;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XPY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(-mp_.t8);
					updateHamiltonian(ham,ix,iy,hopping);

					// yz -> xy, +x-y
					iorb=1; iorb2=2;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XMY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(mp_.t8); // ***********************
					updateHamiltonian(ham,ix,iy,hopping);

					// xy -> yz, +x-y
					iorb=2; iorb2=1;
					ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					isite2=geometry_.getNeighbour(isite,DIR_XMY).first;
					iy=isite2+(iorb2+ispin*mp_.numberOfOrbitals)*volume;
					hopping = signHop(isite)*(-mp_.t8); // ***********************
					updateHamiltonian(ham,ix,iy,hopping);
				}

			/*--------- delta_xy: energy splitting ------------*/
			for (size_t ispin=0;ispin<SPINS;ispin++)
				for (size_t isite = 0; isite < volume; isite++) {
					size_t iorb=2;
					size_t ix=isite+(iorb+ispin*mp_.numberOfOrbitals)*volume;
					ham(ix,ix) += mp_.deltaXY;
				}
		}
	
		RealType hopping(size_t isite,size_t dir2,size_t iorb,size_t iorb2) const
		{
			switch (dir2) {
				case DIR_X:
					if (iorb==0 && iorb2==2) return signHop(isite)*mp_.t7;
					if (iorb==2 && iorb2==0) return signHop(isite)*mp_.t7;
					return 0;
				case DIR_Y:
					if (iorb==1 && iorb2==2) return signHop(isite)*mp_.t7;
					if (iorb==2 && iorb2==1) return signHop(isite)*mp_.t7;
					return 0;
				case DIR_XPY:
					if (iorb==0 && iorb2==2) return signHop(isite)*mp_.t8;
					if (iorb==2 && iorb2==0) return -signHop(isite)*mp_.t8;
					if (iorb==1 && iorb2==2) return signHop(isite)*mp_.t8;
					if (iorb==2 && iorb2==1) return -signHop(isite)*mp_.t8;
					return 0;
				case DIR_XMY:
					if (iorb==0 && iorb2==2) return -signHop(isite)*mp_.t8;
					if (iorb==2 && iorb2==0) return signHop(isite)*mp_.t8;
					if (iorb==1 && iorb2==2) return -signHop(isite)*mp_.t8;
					if (iorb==2 && iorb2==1) return signHop(isite)*mp_.t8;
					return 0;
			}
			return 0;
		}

	private:
		void updateHamiltonian(
				MatrixType& ham,
				size_t ix,
				size_t iy,
				const FieldType& hopping) const
		{
			ham(ix,iy) += hopping;
			ham(iy,ix) += conj(hopping);
		}

		int signHop(size_t i_site) const
		{
			/*------------- sign of t_7 and t_8 ------------*/
			std::vector<int> coor_site(2);
			geometry_.indexToCoor(coor_site,i_site);
			int iSign = coor_site[0] + coor_site[1];
			return (iSign & 1) ? -1 : 1;
//			int sign=1;
//			for(int pos_sign=0;pos_sign<iSign;pos_sign++) sign *= -1;
//			return sign;
		}
/*
		int signHop(size_t isite) const
		{
			return (isite & 1) ? -1 : 1;
		}
*/
		const ParametersModelType& mp_;
		const GeometryType& geometry_;

	}; // class ThreeOrbitalTerms
	
} // namespace Spf

/*@}*/
#endif // THREE_ORBITAL_TERMS_H
