#include "iterator.hpp"

Iterator::Iterator(const Geometry* geom) {
	_geom = geom;
	_value = 0;
	_valid = true;
}

Iterator::Iterator(const Geometry* geom, const index_t& value) {
	_geom = geom;
	_value = value;
	_valid = true;
}

const index_t& Iterator::Value() const {
	return _value;
}

Iterator::operator const index_t&() const {
	return _value;
}

multi_index_t Iterator::Pos() const {
	return _value;
}

void Iterator::First() {
	_value = 0;
}

void Iterator::Next() {
	/*multi_index_t size = _geom->Size();
	if ( _value + 1 >= size[0]*size[1]-1)
		_valid = false;
	else
		_value += 1;*/
}

bool Iterator::Valid() const {
	return _valid;
}

Iterator Iterator::Left() const {
/*	if ( (_value - 1) % _geom->Size()[0] == 0)
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value - 1);*/
}

Iterator Iterator::Right() const {
	/*if ( _value % _geom->Size()[0] == 0)
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value + 1);*/
}

Iterator Iterator::Top() const {
	//return Iterator(_geom, _value - _geom->Size()[0]);
}

Iterator Iterator::Down() const {
	//return Iterator(_geom, _value + _geom->Size()[0]);
}

InteriorIterator::InteriorIterator(const Geometry* geom) : Iterator(geom) {

}

void InteriorIterator::First() {
}

void InteriorIterator::Next() {
}

BoundaryIterator::BoundaryIterator(const Geometry* geom) : Iterator(geom) {
}

void BoundaryIterator::SetBoundary(const index_t& boundary) {
	_boundary = boundary;
}

void BoundaryIterator::First() {
}

void BoundaryIterator::Next() {
}
