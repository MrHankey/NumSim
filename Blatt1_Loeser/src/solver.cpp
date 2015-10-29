#include "solver.hpp"

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
	norm  = 0.5*(dx^2+dy^2);

	real_t pDiff = piij+pijj-4*pij+pj+pi;

	return 1/norm * (rhs-pDiff);
}

SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	_omega = omega;
}

SOR::~SOR() {
}

real_t SOR::Cycle(Grid* grid, const Grid* rhs) const {
	Iterator it;

	real_t n, res, dx, dy, norm;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*dx^2*dy^2/(dx^2+dy^2);

	while(it.Valid()) {
		n++;

		real_t pij, pi, pj, pijj, piij, A, B, C, corr;

		pij  = grid->Cell(it);
		pi   = grid->Cell(it.Down());
		pijj = grid->Cell(it.Top());
		pj   = grid->Cell(it.Left());
		piij = grid->Cell(it.Right());

		A    = (pj+pj)/dx^2;
		B    = (pi+pijj)/dy^2;
		C    = 1/norm*pij;

		corr = A+B-C-rhs;
		res  = 0;

		pij  = pij + _omega * norm * corr;
		res += localRes(it,grid,rhs);

		it.Next();

	}
	return res/n;
}
