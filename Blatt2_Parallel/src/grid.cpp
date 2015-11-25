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

#include "grid.hpp"

#include <iostream>
#include <cmath>

#include "iterator.hpp"
#include "geometry.hpp"

using namespace std;

/// Constructs a grid based on a geometry
//  @param geom    Geometry information
Grid::Grid(const Geometry* geom) {
	multi_real_t offset;
	offset[0] = 0;
	offset[1] = 0;

	_data     = new real_t[geom->Size()[0]*geom->Size()[1]];
	_geom     = geom;
	_offset   = offset;
}

/// Constructs a grid based on a geometry with an offset
//  @param geom    Geometry information
//  @param offset  distance of staggered grid point to cell's anchor point
Grid::Grid(const Geometry* geom, const multi_real_t& offset) {
	_data = new real_t[geom->Size()[0]*geom->Size()[1]];
	_offset = offset;
	_geom = geom;
}

/// Destructor: deletes the grid
Grid::~Grid() {delete[] _data;}

/// Initializes the grid with a given value
//  @param value  value to initialize grid
void Grid::Initialize(const real_t& value) {
	for (index_t i = 0; i< _geom->Size()[0]*_geom->Size()[1];i++){
		_data[i] = value;
	}
}

/// Write access to the grid cell at position [it]
//  @param it  iterator position
real_t& Grid::Cell(const Iterator& it) {
	return _data[it];
}

/// Read access to the grid cell at position [it]
//  @param it  iterator position
const real_t& Grid::Cell(const Iterator& it) const {
	return _data[it];
}

/// Returns a bilinear interpolation at an arbitrary position
//  @param  pos  Position given between 0 to 1
real_t Grid::Interpolate(const multi_real_t& pos) const {
	// Initialize
	real_t xP, x1, x2, yP, y1, y2;
	real_t q11, q12, q21, q22;
	real_t r1, r2;

	// Get cell where position is at - row, column and cell id
	index_t numCol = (index_t)((pos[0] - _offset[0])*(_geom->Size()[0] - 2) + 1);
	index_t numRow = (index_t)((pos[1] - _offset[1])*(_geom->Size()[1] - 2) + 1);
	index_t cellID = (index_t)((numRow)*_geom->Size()[0])+numCol;

	// Initialize iteraytor at starting position of the cell id
	Iterator it = Iterator(_geom,cellID);

	// Get the four surrounding cell values
	q11 = _data[it];
	q12 = _data[it.Right()];
	q21 = _data[it.Top()];
	q22 = _data[it.Top().Right()];

	// Get x positions of the surrounding cells
	xP = pos[0];
	x1 = (real_t)(it.Pos()[0])/(real_t)(_geom->Size()[0] - 1);
	x2 = x1 + _geom->Mesh()[0];

	// Get y positions of the surrounding cells
	yP = pos[1];
	y1 = (real_t)(it.Pos()[1])/(real_t)(_geom->Size()[1] - 1);
	y2 = y1 + _geom->Mesh()[1];

	// Interpolate in x direction for each two cells
	r1 = (x2-xP)/(x2-x1)*q11 + (xP-x1)/(x2-x1)*q12;
	r2 = (x2-xP)/(x2-x1)*q21 + (xP-x1)/(x2-x1)*q22;

	// Interpolate in y direction with x interpolation
	return (y2-yP)/(y2-y1)*r1 + (yP-y1)/(y2-y1)*r2;
}

/* Calculate differences */
/// Computes the left-sided difference quatient in x-dim at [it]
//  @param it  iterator position
real_t Grid::dx_l(const Iterator& it) const {
	return (_data[it]-_data[it.Left()])/_geom->Mesh()[0];
}

/// Computes the right-sided difference quatient in x-dim at [it]
//  @param it  iterator position
real_t Grid::dx_r(const Iterator& it) const {
	return (_data[it.Right()]-_data[it])/_geom->Mesh()[0];
}

/// Computes the left-sided difference quatient in y-dim at [it]
//  @param it  iterator position
real_t Grid::dy_l(const Iterator& it) const {
	return (_data[it]-_data[it.Down()])/_geom->Mesh()[1];
}

/// Computes the right-sided difference quatient in x-dim at [it]
//  @param it  iterator position
real_t Grid::dy_r(const Iterator& it) const {
	return (_data[it.Top()]-_data[it])/_geom->Mesh()[1];
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
//  @param it  iterator position
real_t Grid::dxx(const Iterator& it) const {
	return (_data[it.Right()]-2*_data[it]+_data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
//  @param it  iterator position
real_t Grid::dyy(const Iterator& it) const {
	return (_data[it.Top()]-2*_data[it]+_data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
}

/* Donor cell method */
/// Computes u*du/dx with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const {
	real_t A = (_data[it]+_data[it.Right()])/2;
	real_t B = (_data[it.Left()]+_data[it])/2;
	real_t C = (_data[it]-_data[it.Right()])/2;
	real_t D = (_data[it.Left()]-_data[it])/2;
	return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[0];
}

/// Computes v*du/dy with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
//  @param v      get grid for velocity v
real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha,const Grid* v) const {
	real_t A = (v->Cell(it)+v->Cell(it.Right()))/2;
	real_t B = (_data[it]+_data[it.Top()])/2;
	//std::cout<<"B: "<<B;
	real_t C = (v->Cell(it.Down())+ v->Cell(it.Down().Right()))/2;
	real_t D = (_data[it]+_data[it.Down()] )/2;
	real_t E = (_data[it]-_data[it.Top()])/2;
	real_t F = (_data[it.Down()]-_data[it])/2;
	return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[1];
}

/// Computes u*dv/dx with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
//  @param u      get grid for velocity u
real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha,const Grid* u) const {
	real_t A = (u->Cell(it)+u->Cell(it.Top()))/2;
	real_t B = (_data[it]+_data[it.Right()])/2;
	real_t C = (u->Cell(it.Left())+ u->Cell(it.Left().Top()))/2;
	real_t D = (_data[it]+_data[it.Left()] )/2;
	real_t E = (_data[it]-_data[it.Right()])/2;
	real_t F = (_data[it.Left()]-_data[it] )/2;
	return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[0];

}

/// Computes v*dv/dy with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
	real_t A = (_data[it]+_data[it.Top()])/2;
	real_t B = (_data[it.Down()]+_data[it])/2;
	real_t C = (_data[it]-_data[it.Top()])/2;
	real_t D = (_data[it.Down()]-_data[it])/2;
	return ( (A*A - B*B) + alpha*( fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[1];
}


/// Returns the maximal value of the grid
real_t Grid::Max() const {
	real_t max = _data[0];
	for (index_t i = 1; i< _geom->Size()[0]*_geom->Size()[1];i++){
		if (_data[i] >= max)
		max = _data[i];
	}
	return max;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
	real_t min = _data[0];
	for (index_t i = 1; i< _geom->Size()[0]*_geom->Size()[1];i++){
		if (_data[i] <= min)
		min = _data[i];
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
	real_t max = this->Max();
	real_t min = this->Min();
	if ((max+min) > 0)
		return max;
	else
		return min;
}

/// Returns a pointer to the raw data
real_t* Grid::Data() {
	return _data;
}

const multi_real_t& Grid::getOffset() const {
	return _offset;
}

const Geometry* Grid::getGeometry() const {
	return _geom;
}
