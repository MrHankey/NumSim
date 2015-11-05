#include "solver.hpp"
#include "grid.hpp"
#include "geometry.hpp"

#include <cmath>

Solver::Solver(const Geometry* geom) {
	_geom = geom;
}

Solver::~Solver() {
}

real_t Solver::localRes(const Iterator& it, const Grid* grid,
		const Grid* rhs) const {
	/*real_t dx, dy, pij, pij_d, pij_l, pij_t, pij_r, norm, pDiff;

	dx    = _geom->Mesh()[0];
	dy    = _geom->Mesh()[1];

	pij   = grid->Cell(it);
	pij_d = grid->Cell(it.Down());
	pij_t = grid->Cell(it.Top());
	pij_l = grid->Cell(it.Left());
	pij_r = grid->Cell(it.Right());
	norm  = 0.5*(dx*dx+dy*dy);

	pDiff = (pij_l - 2.0*pij + pij_r)/(dx*dx) + (pij_d - 2.0*pij + pij_t)/(dy*dy);

	return (rhs->Cell(it)-pDiff);*/
	return abs(grid->dxx(it) + grid->dyy(it) - rhs->Cell(it));
}

SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	//include both h_x and h_y in calculation
	_omega = omega;//2.0/(1.0 + std::sin(M_PI*geom->Mesh()[0]));//omega;
}

SOR::~SOR() {
}

real_t SOR::Cycle(Grid* grid, const Grid* rhs) const {
	InteriorIterator it = InteriorIterator(_geom);

	real_t n, res, dx, dy, norm;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	while(it.Valid()) {
		n++;

		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		pij  	= grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t 	= grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r 	= grid->Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		grid->Cell(it)  = (1-_omega)*pij + _omega * norm * corr;
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		it.Next();

	}
	return res/n;
}
