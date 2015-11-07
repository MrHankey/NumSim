#include "solver.hpp"
#include "grid.hpp"
#include "geometry.hpp"

#include <cmath>

/// Class Constuctor
// @param geom  get field geometry
Solver::Solver(const Geometry* geom) {
	_geom = geom;
}

Solver::~Solver() {
}

/// Calculate Residual
// @param  it        calls iterator
// @param  grid      get grid values
// @return localRes  return step residual
real_t Solver::localRes(const Iterator& it, const Grid* grid,
		const Grid* rhs) const {
	return fabs(grid->dxx(it) + grid->dyy(it) - rhs->Cell(it));
}

/// SOR Constructor
// @param geom   get geometry
// @param omega  get scaling factor
SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	_omega = omega;
}

SOR::~SOR() {
}

/// Implement SOR method
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t SOR::Cycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy, norm;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	// Cycle through all cells
	while(it.Valid()) {
		// Initialize
		n++;
		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		// define variables
		pij  	= grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t 	= grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r 	= grid->Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		// Save new p_ij
		grid->Cell(it)  = (1-_omega)*pij + _omega * norm * corr;
		// Calculate and summ residual
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		// Next cell
		it.Next();
	}

	// Norm residual
	return res/n;
}
