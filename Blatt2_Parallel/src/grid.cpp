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

  if (true) {
	  //TODO do proper interpolation
	  index_t cell_x = (index_t)ceil((pos[0]/_geom->Mesh()[0]));
	  index_t cell_y = (index_t)ceil((pos[1]/_geom->Mesh()[1]));

	  return _data[cell_y*_geom->Size()[0] + cell_x];
  } else {
	//subtracts offsets
	multi_real_t p_off;
	p_off[0] = pos[0] - _offset[0];
	p_off[1] = pos[1] - _offset[1];

	//relative position of x and y in physical grid
	real_t relx = p_off[0]/_geom->Length()[0];
	real_t rely = p_off[1]/_geom->Length()[1];

	//calculates x and y coordinate of cell
	index_t cell1 = floor(relx*(_geom->Size()[0] - 2.0)+1);
	index_t cell2 = floor(rely*(_geom->Size()[1] - 2.0)+1);

	//transfers x and y coordinates of cell to iterator
	index_t value = cell2 * _geom->Size()[0] + cell1;
	Iterator it(_geom, value);

	//x and y values for interpolation
	real_t x1 = ((real_t)it.Pos()[0] -1)* _geom->Mesh()[0] + _offset[0];
	real_t x2 = x1 + _geom->Mesh()[0];
	real_t y1 = ((real_t)it.Pos()[1]-1)* _geom->Mesh()[1] + _offset[1];
	real_t y2 = y1 + _geom->Mesh()[0];

	real_t data = _data[it];
	real_t data_t = _data[it.Top()];
	real_t data_r = _data[it.Right()];
	real_t data_rt = _data[it.Right().Top()];

	if ( _geom->_b->Cell(it) == 1 ) data = 0;
	if ( _geom->_b->Cell(it.Top()) == 1 ) data_t = 0;
	if ( _geom->_b->Cell(it.Right()) == 1 ) data_r = 0;
	if ( _geom->_b->Cell(it.Right().Top()) == 1 ) data_rt = 0;


	//interpolate in x direction
	real_t r1 = data       * ((x2 - pos[0])/_geom->Mesh()[0]) + data_r * ((pos[0] - x1)/_geom->Mesh()[0]);
	real_t r2 = data_t * ((x2 - pos[0])/_geom->Mesh()[0]) + data_rt * ((pos[0] - x1)/_geom->Mesh()[0]);

	//interpolate in y direction
	return r1 * ((y2 - p_off[1])/_geom->Mesh()[1]) + r2 * ((p_off[1] - y1)/_geom->Mesh()[1]);
  }
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
	real_t A = (_data[it]+_data[it.Right()])/2.0;
	real_t B = (_data[it.Left()]+_data[it])/2.0;
	real_t C = (_data[it]-_data[it.Right()])/2.0;
	real_t D = (_data[it.Left()]-_data[it])/2.0;
	return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[0];
}

/// Computes v*du/dy with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
//  @param v      get grid for velocity v
real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha,const Grid* v) const {
	real_t A = (v->Cell(it)+v->Cell(it.Right()))/2.0;
	real_t B = (_data[it]+_data[it.Top()])/2.0;
	//std::cout<<"B: "<<B;
	real_t C = (v->Cell(it.Down())+ v->Cell(it.Down().Right()))/2.0;
	real_t D = (_data[it]+_data[it.Down()] )/2.0;
	real_t E = (_data[it]-_data[it.Top()])/2.0;
	real_t F = (_data[it.Down()]-_data[it])/2.0;
	return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[1];
}

/// Computes u*dv/dx with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
//  @param u      get grid for velocity u
real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha,const Grid* u) const {
	real_t A = (u->Cell(it)+u->Cell(it.Top()))/2.0;
	real_t B = (_data[it]+_data[it.Right()])/2.0;
	real_t C = (u->Cell(it.Left())+ u->Cell(it.Left().Top()))/2.0;
	real_t D = (_data[it]+_data[it.Left()] )/2.0;
	real_t E = (_data[it]-_data[it.Right()])/2.0;
	real_t F = (_data[it.Left()]-_data[it] )/2.0;
	return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[0];

}

/// Computes v*dv/dy with the donor cell method
//  @param it     iterator position
//  @param alpha  parameter for donor cell method
real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
	real_t A = (_data[it]+_data[it.Top()])/2.0;
	real_t B = (_data[it.Down()]+_data[it])/2.0;
	real_t C = (_data[it]-_data[it.Top()])/2.0;
	real_t D = (_data[it.Down()]-_data[it])/2.0;
	return ( (A*A - B*B) + alpha*( fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[1];
}


/// Returns the maximal value of the grid
real_t Grid::Max() const {
	Iterator it = Iterator(_geom);

	real_t max = Cell(it);
	while (it.Valid())
	{
		if (Cell(it) >= max)
			max = Cell(it);

		it.Next();
	}
	return max;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
	Iterator it = Iterator(_geom);

	real_t min = Cell(it);
	while (it.Valid())
	{
		if (Cell(it) <= min)
			min = Cell(it);

		it.Next();
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
	real_t max = this->Max();
	real_t min = this->Min();
	if ((max+min) > 0)
		return fabs(max);
	else
		return fabs(min);
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
