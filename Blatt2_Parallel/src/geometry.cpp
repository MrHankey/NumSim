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

#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;

/// Constructs a default geometry:
/// driven cavity with 128 x 128 grid, no-slip boundary conditions
Geometry::Geometry() {
	// Set starting velocities
	_velocity[0] = 1;
	_velocity[1] = 0;
	_pressure    = 0;


	// Number of cells in one line
	_size[0] = 128;
	_size[1] = 128;

	// Length of driven cavity
	_length[0] = 1;
	_length[1] = 1;

	// Print vars
	cout << "Loaded default geometry definition." << endl;
}

/// Loads a geometry from a file
//  @param file  txt file to load data from
void Geometry::Load(const char* file) {
	// Initialize
	ifstream fileStream(file);

	// Cath loading error
	if ( fileStream.fail() ) {
		// Failed to load
		cout << "Couldn't load geometry from " << file << endl;
		return;
	} else {
		//Load values
		fileStream >> _velocity[0];
		fileStream >> _velocity[1];
		fileStream >> _pressure;
		fileStream >> _size[0];
		fileStream >> _size[1];
		fileStream >> _length[0];
		fileStream >> _length[1];

		// Success
		cout << "Loaded geometry definitions from " << file << "." << endl;
	}
	Initialize();
}

void Geometry::Initialize()
{
	// Cell Length
	_h[0] = _length[0]/(_size[0]);
	_h[1] = _length[1]/(_size[1]);

	_size[0] += 2;
	_size[1] += 2;

	printf(" local_siz: %i \n", _size[0]);
}

/// Prints Parameters
void Geometry::PrintVariables(){
	cout << "vel_x: "    << _velocity[0]  << endl;
	cout << "vel_y: "    << _velocity[1]  << endl;
	cout << "p: "        << _pressure     << endl;
	cout << "size_x: "   << _size[0]      << endl;
	cout << "size_y: "   << _size[1]      << endl;
	cout << "length_x: " << _length[0]    << endl;
	cout << "length_y "  << _length[1]    << endl;
	cout << "h_x: "      << _h[0]         << endl;
	cout << "h_y: "      << _h[1] << endl << endl;

}

/* Getter functions */
/// Returns the number of cells in each dimension
const multi_index_t& Geometry::Size()   const {
	return _size;
}
/// Returns the length of the domain
const multi_real_t&  Geometry::Length() const
{
	return _length;
}
/// Returns the meshwidth
const multi_real_t&  Geometry::Mesh()   const {return _h;}


/// Updates the boundary velocity field u
//  @param u  get grid of u
void Geometry::Update_U(Grid* u) const {
	// Initialize
	real_t velocity = _velocity[0];
	BoundaryIterator it = BoundaryIterator(this);


		it.SetBoundary(it.boundaryBottom);
		while ( it.Valid() ) {
			u->Cell(it) = -1*u->Cell(it.Top());
			it.Next();
		}



		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			u->Cell(it) = 0;
			u->Cell(it.Left()) = 0; // hier evtl deleten
			it.Next();
		}


		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			u->Cell(it) = 2*velocity - u->Cell(it.Down()) ;
			it.Next();
		}

		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			u->Cell(it) = 0;
			it.Next();
		}

}

/// Updates the boundary velocity field v
//  @param v  get grid of v
void Geometry::Update_V(Grid* v) const {
	// Initialize
	BoundaryIterator it = BoundaryIterator(this);


		it.SetBoundary(it.boundaryBottom);
		while ( it.Valid() ) {
			v->Cell(it) = 0;
			it.Next();
		}

		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			v->Cell(it) = -1*v->Cell(it.Left());
			it.Next();
		}

		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			v->Cell(it.Down()) = 0; // hier evtl deleten
			v->Cell(it) = 0;
			it.Next();
		}

		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			v->Cell(it) = -1*v->Cell(it.Right());
			it.Next();
		}

}

/// Updates the boundary velocity field p
//  @param p  get grid of p
void Geometry::Update_P(Grid* p) const {
	// Initialize
	BoundaryIterator it = BoundaryIterator(this);
	it.SetBoundary(it.boundaryBottom);


		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Top());
			it.Next();
		}

		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Left());
			it.Next();
		}

		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Down());
			it.Next();
		}


		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Right());
			it.Next();
		}

}


