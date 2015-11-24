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

#include "iterator.hpp"

#include "geometry.hpp"

#include <iostream>
/// Constructs a new Iterator depending on a geometry
//  @param geom  get geometry
Iterator::Iterator(const Geometry* geom) {
	_geom = geom;
	First();
}

/// Constructs a new Iterator on a geometry with a defined starting value
//  @param geom   get geometry
//  @param value  starting value for iteration
Iterator::Iterator(const Geometry* geom, const index_t& value) {
	_geom = geom;
	_value = value;
	_valid = true;
}

/// Returns the current position value
const index_t& Iterator::Value() const {
	return _value;
}

/// Cast operator to convert Iterators to integers
Iterator::operator const index_t&() const {
	return _value;
}

/// Returns the position coordinates
multi_index_t Iterator::Pos() const {
	multi_index_t pos;
	pos[0] = _value % _geom->Size()[0];
	pos[1] = ( _value - pos[0] ) / _geom->Size()[0];
	return pos;
}

/// Sets the iterator to the first element
void Iterator::First() {
	_value = 0;
	_valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
	multi_index_t size = _geom->Size();
	if ( _value + 1 > size[0]*size[1]-1)
		_valid = false;
	else
		_value += 1;
}

/// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
	return _valid;
}

/// Returns an Iterator that is located left from this one.
//  On left boundary iterator returns the same cell
Iterator Iterator::Left() const {
	if ( _value % _geom->Size()[0] == 0 )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value - 1);
}

/// Returns an Iterator that is located right from this one
//  On right boundary iterator returns the same cell
Iterator Iterator::Right() const {
	if ( (_value + 1) % (_geom->Size()[0]) == 0)
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value + 1);
}

/// Returns an Iterator that is located above this one
//  On top boundary iterator returns the same cell
Iterator Iterator::Top() const {
	index_t newIndex = _value + _geom->Size()[0];
	if ( newIndex > (_geom->Size()[0]*_geom->Size()[1] - 1) )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, newIndex);
}

/// Returns an Iterator that is located below this one
//  On bottom boundary iterator returns the same cell
Iterator Iterator::Down() const {
	int newIndex = _value - _geom->Size()[0];
	if ( newIndex < 0 )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, newIndex);
}


/* Interior iterator */
/// Construct a new InteriorIterator
//  @param geom   get geometry
InteriorIterator::InteriorIterator(const Geometry* geom) : Iterator(geom) {
	First();
}

/// Sets the iterator to the first element
void InteriorIterator::First() {
	_value = _geom->Size()[0] + 1;
	_valid = true;
	//std::cout <<_value<<" ";
}

/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
	if ( (_value + 2) % (2*_geom->Size()[0]) == 0 ){
		_value += _geom->Size()[0];
		if ( (_value) >= (_geom->Size()[0]*(_geom->Size()[1] - 1)))
			_valid = false;
	}
	else if ( (_value - 1) % (2*_geom->Size()[0]) == 0 ){
		_value += _geom->Size()[0];
		if ( (_value) >= (_geom->Size()[0]*(_geom->Size()[1] - 1)))
			_valid = false;
	}
	else if (((index_t)(_value/_geom->Size()[0])) % 2 == 1)
		_value++;
	else if (((index_t)(_value/_geom->Size()[0])) % 2 == 0)
		_value--;
	else
		std::cout<<_value<<" erzeugt Fehler!";

}


/* Boundary interator */
/// Constructs a new BoundaryIterator
//  @param geom   get geometry
BoundaryIterator::BoundaryIterator(const Geometry* geom) : Iterator(geom) {
	_boundary = 0;
	First();
}

/// Sets the boundary to iterate
//  @param boundary  get which boundary to iterate
//  0 - Bottom boundary
//  1 - Right  boundary
//  2 - Top    boundary
//  3 - Left   boundary
void BoundaryIterator::SetBoundary(const index_t& boundary) {
	_boundary = boundary;
	First();
}

/// Sets the iterator to the first element
void BoundaryIterator::First() {
	_valid = true;

	// Bottom boundary
	if (_boundary == 0 ) {
		_value = 0;
	}
	// Right boundary
	else if ( _boundary == 1 ) {
		_value = _geom->Size()[0] - 1;
	}
	// Top boundary
	else if ( _boundary == 2 ) {
		_value = _geom->Size()[0]*(_geom->Size()[1]-1);
	}
	// Left boundary
	else if ( _boundary == 3 ) {
		_value = 0;
	}
}

/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
	// Bottom boundary
	if (_boundary == 0 ) {
		_value++;
		if ( _value >= _geom->Size()[0] ) {
			_valid = false;
		}
	}
	// Right boundary
	else if ( _boundary == 1 ) {
		_value += _geom->Size()[0];
		if ( _value >= _geom->Size()[0]*(_geom->Size()[1]) ) {
			_valid = false;
		}
	}
	// Top boundary
	else if ( _boundary == 2 ) {
		_value++;
		if ( _value >= _geom->Size()[0]*(_geom->Size()[1]) ) {
			_valid = false;
		}
	}
	// Left boundary
	else if ( _boundary == 3 ) {
		_value += _geom->Size()[0];
		if ( _value >= _geom->Size()[0]*(_geom->Size()[1]) ) {
			_valid = false;
		}
	}
}
