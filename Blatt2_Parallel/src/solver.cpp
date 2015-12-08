/*
 * Copyright (C) 2015  Raphael Leiteriz, Sebastian Reuschen, Hamzeh Kraus
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "solver.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

#include <cmath>
#include <iostream>

/// Constructor of the abstract Solver class
//  @param geom  get field geometry
Solver::Solver(const Geometry* geom) {
	_geom = geom;
}

/// Destructor of the Solver Class
Solver::~Solver() {}

/// Returns the residual at [it] for the pressure-Poisson equation
//  @param  it        calls iterator
//  @param  grid      get grid values
//  @return localRes  return residual
real_t Solver::localRes(const Iterator& it, const Grid* grid,
		const Grid* rhs) const {
	return fabs(grid->dxx(it) + grid->dyy(it) - rhs->Cell(it));
}

/* SOR solver */
/// Concrete SOR solver
//  @param geom   get geometry
//  @param omega  get scaling factor
SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	_omega = omega;
}

/// Destructor
SOR::~SOR() {}


/// Returns the total residual and executes a solver cycle of SOR
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

//---------------------------------------------------------------------------------------------------

/* SOR solver */
/// Concrete SOR solver
//  @param geom   get geometry
//  @param omega  get scaling factor
RedOrBlackSOR::RedOrBlackSOR(const Geometry* geom, const real_t& omega) : SOR(geom,omega) {
	_geom  = geom;
	_omega = omega;
}

/// Destructor
RedOrBlackSOR::~RedOrBlackSOR() {}

/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t RedOrBlackSOR::RedCycle(Grid* grid, const Grid* rhs) const {
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
		if(res<=lRes){
			res=lRes;//max statt bisher mean
		}
		//res += lRes;

		// Next cell
		it.Next();
		it.Next();
	}

	// Norm residual
	//return res/n;
	return res;
}


/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t RedOrBlackSOR::BlackCycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy, norm;
	it.Next();
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

		if(res<lRes){
			res=lRes;//max statt bisher mean
		}
		//res += lRes;

		// Next cell
		it.Next();
		it.Next();
	}

	// Norm residual
	//return res/n; //max statt bisher mean
	return res;
}


