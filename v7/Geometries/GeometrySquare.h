
/** \ingroup SPF */
/*@{*/

/*! \file GeometrySquare.h
 *
 * A square lattice
 *
 */
#ifndef GEOM_SQUARE_H
#define GEOM_SQUARE_H
#include <utility>
#include "String.h"
#include <vector>
#include "Matrix.h" // in PsimagLite
#include "ProgramGlobals.h"

namespace Spf {
template<typename FieldType_>
class GeometrySquare {

public:

	typedef std::pair<SizeType,SizeType> PairType;

	GeometrySquare(SizeType l) : l_(l),volume_(l*l)
	{
		buildNeighbors();
	}

	SizeType distances() const { return neighbors_.size(); }

	SizeType z(SizeType distance=1) const
	{
		return neighbors_[distance-1].n_col();
	}

	// j-th neighbor of i at distance (starts from 1 for compatibility)
	PairType neighbor(SizeType i,SizeType j,SizeType distance=1) const
	{
		return neighbors_[distance-1](i,j);
	}

	PairType getNeighbour(SizeType i,SizeType dir) const
	{
		SizeType j = dir*2;
		SizeType distance = 1;
		if (j>=4) {
			distance++;
			j -= 4;
		}
		return neighbors_[distance-1](i,j);
	}

	SizeType volume() const { return volume_; }

	SizeType add(SizeType ind,SizeType ind2) const
	{
		PsimagLite::Vector<int>::Type x(2),y(2);
		indexToCoor(x,ind);
		indexToCoor(y,ind2);
		for (SizeType i=0;i<x.size();i++) {
			x[i] += y[i];
			g_pbc(x[i],l_);
		}
		return g_index(x);
	}

	SizeType dim() const { return 2; }

	SizeType length() const { return l_; }

	PsimagLite::String name() const { return "square"; }

	template<typename SomeVectorType>
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,void>::Type
	indexToCoor(SomeVectorType& v,SizeType i) const
	{
		SizeType lx = l_;
		v[0] = i%lx;
		v[1] = SizeType(i/lx);
	}

	SizeType coorToIndex(SizeType x,SizeType y) const
	{
		SizeType lx = l_;
		return x + y*lx;
	}

	// direction connecting i and j
	// assume that i and j are n-neighbors or next n-neighbors
	int getDirection(SizeType ind,SizeType jnd) const
	{
		PsimagLite::Vector<int>::Type vi(2),vj(2);
		indexToCoor(vi,ind);
		indexToCoor(vj,jnd);

		for (SizeType i=0;i<vi.size();i++) {
			vi[i] -= vj[i];
		}

		if (vi[0] == 0 && periodicEqualTo(vi[1],1)) return ProgramGlobals::DIRY;
		if (vi[0] == 0 && periodicEqualTo(vi[1],-1)) return ProgramGlobals::DIRY;

		if (vi[1] == 0 &&  periodicEqualTo(vi[0],1)) return ProgramGlobals::DIRX;
		if (vi[1] == 0 &&  periodicEqualTo(vi[0],-1)) return ProgramGlobals::DIRX;

		if (periodicEqualTo(vi[0],1) && periodicEqualTo(vi[1],1)) return ProgramGlobals::DIRXPY;
		if (periodicEqualTo(vi[0],-1) && periodicEqualTo(vi[1],-1)) return ProgramGlobals::DIRXPY;

		if (periodicEqualTo(vi[0],1) && periodicEqualTo(vi[1],-1)) return ProgramGlobals::DIRXMY;
		if (periodicEqualTo(vi[0],-1) && periodicEqualTo(vi[1],1)) return ProgramGlobals::DIRXMY;

		return -1;
	}

	SizeType scalarDirection(SizeType,SizeType) const
	{
		throw PsimagLite::RuntimeError("scalarDirection unimplemented\n");
	}

	template<typename SomeVectorType>
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,SizeType>::Type
	coor2Index(const SomeVectorType& v) const
	{
		throw PsimagLite::RuntimeError("coor2Index unimplemented\n");
	}

	template<typename T>
	friend std::ostream& operator<<(std::ostream& os,
	                                const GeometrySquare<T>& g);

private:

	bool periodicEqualTo(int a,int b) const
	{
		if (a==b) return true;
		int c = -b*(l_-1);
		if (a==c) return true;
		return false;
	}

	void buildNeighbors()
	{
		neighborsAt1();
		neighborsAt2();
	}

	void neighborsAt1()
	{
		int lx = l_;
		int ly = l_;
		int zz = 0;
		//PairType zeroVal(0,0);
		PsimagLite::Matrix<PairType> matrix(lx*ly,4);
		for (int y=0;y<ly;y++) {
			for (int x=0;x<lx;x++) {
				SizeType i = x + y*lx;
				int xx=x+1;
				int yy=y;
				SizeType counter = 0;
				/*if (g_pbc(xx,lvector[0])) {
					border.push_back(0);
				} else {
					border.push_back(-1);
				}*/
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRX);
				xx=x-1;
				/*if (g_pbc(xx,lvector[0])) {
					border.push_back(0);
				} else {
					border.push_back(-1);
				}*/
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRX);
				xx=x; yy=y+1;
				/*if (g_pbc(yy,lvector[1])) {
					border.push_back(1);
				} else {
					border.push_back(-1);
				}*/
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRY);
				yy=y-1;
				/*if (g_pbc(yy,lvector[1])) {
					border.push_back(1);
				} else {
					border.push_back(-1);
				}*/
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRY);
			}
		}
		neighbors_.push_back(matrix);
	}

	void neighborsAt2()
	{
		int lx = l_;
		int ly = l_;
		int zz = 0;
		//PairType zeroVal(0,0);
		PsimagLite::Matrix<PairType> matrix(lx*ly,4);
		for (int y=0;y<ly;y++) {
			for (int x=0;x<lx;x++) {
				SizeType i = x + y*lx;
				int xx=x+1;
				int yy=y+1;
				SizeType counter=0;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRXPY);
				xx=x-1; yy=y-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRXPY);
				xx=x-1; yy=y+1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRXMY);
				xx=x+1; yy=y-1;
				matrix(i,counter++) = PairType(g_index(xx,yy,zz),ProgramGlobals::DIRXMY);
			}
		}
		neighbors_.push_back(matrix);
	}

	bool g_pbc(int& x, SizeType l) const
	{
		int L = l;
		bool r=false;
		if (x<0) r=true;
		if (x>=L) r=true;
		while (x<0) x+=L;
		while (x>=L) x-=L;
		return r;
	}

	SizeType g_index(PsimagLite::Vector<int>::Type& x) const
	{
		int zz=0;
		return g_index(x[0],x[1],zz);
	}

	SizeType g_index(int& x,int& y,int&) const
	{
		SizeType lx = l_;
		SizeType ly = l_;
		g_pbc(x,lx);
		g_pbc(y,ly);
		//g_pbc(z,lz);
		return x+y*lx; //+z*L*L;
	}

	SizeType l_;
	SizeType volume_;
	PsimagLite::Vector<PsimagLite::Matrix<PairType> >::Type neighbors_;
};

template<typename T>
std::ostream& operator<<(std::ostream& os,const GeometrySquare<T>& g)
{
	for (SizeType i=0;i<g.neighbors_.size();i++) {
		os<<"#i="<<i<<"\n";
		for (SizeType k=0;k<g.neighbors_[i].n_row();k++) {
			for (SizeType l=0;l<g.neighbors_[i].n_col();l++) {
				os<<g.neighbors_[i](k,l).first<<" ";
			}
			os<<"\n";
		}
	}
	return os;
}
} // namespace Spf

/*@}*/
#endif

