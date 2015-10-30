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
	real_t dx, dy, pij, pi, pj, pijj, piij, norm, pDiff;

	dx    = _geom->Mesh()[0];
	dy    = _geom->Mesh()[1];

	pij   = grid->Cell(it);
	pi    = grid->Cell(it.Down());
	pijj  = grid->Cell(it.Top());
	pj    = grid->Cell(it.Left());
	piij  = grid->Cell(it.Right());
	norm  = 0.5*(dx*dx+dy*dy);

	pDiff = piij+pijj-4*pij+pj+pi;

	return (1/norm * (rhs->Cell(it)-pDiff));
}

SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	//include both h_x and h_y in calculation
	_omega = 2.0/(1.0 + std::sin(M_PI*geom->Mesh()[0]));//omega;
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

		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, C, corr;

		pij  = grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t = grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r = grid->Cell(it.Right());

		A    = (pij_l+pij_r)/dx*dx;
		B    = (pij_d+pij_t)/dy*dy;
		C    = 1/norm*pij;

		corr = norm * ( A+B-C-rhs->Cell(it) );

		grid->Cell(it)  = pij + _omega * corr;
		real_t lRes = localRes(it,grid,rhs);
		res += lRes*lRes;

		it.Next();

	}
	return res/n;
}
