#include "iterator.hpp"

Iterator::Iterator(const Geometry* geom) {
	_geom = geom;
	First();
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
	multi_index_t pos;
	pos[0] = _value % _geom->Size()[0];
	pos[1] = ( _value - pos[0] ) / _geom->Size()[0];
	return pos;
}

void Iterator::First() {
	_value = 0;
	_valid = true;
}

void Iterator::Next() {
	multi_index_t size = _geom->Size();
	if ( _value + 1 > size[0]*size[1]-1)
		_valid = false;
	else
		_value += 1;
}

bool Iterator::Valid() const {
	return _valid;
}

Iterator Iterator::Left() const {
	if ( _value % _geom->Size()[0] == 0 )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value - 1);
}

Iterator Iterator::Right() const {
	if ( (_value + 1) % (_geom->Size()[0]) == 0)
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, _value + 1);
}

Iterator Iterator::Top() const {
	index_t newIndex = _value + _geom->Size()[0];
	if ( newIndex > (_geom->Size()[0]*_geom->Size()[1] - 1) )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, newIndex);
}

Iterator Iterator::Down() const {
	int newIndex = _value - _geom->Size()[0];
	if ( newIndex < 0 )
		return Iterator(_geom, _value);
	else
		return Iterator(_geom, newIndex);
}

InteriorIterator::InteriorIterator(const Geometry* geom) : Iterator(geom) {
	First();
}

void InteriorIterator::First() {
	_value = _geom->Size()[0] + 1;
	_valid = true;
	//std::cout <<_value<<" ";
}

void InteriorIterator::Next() {
	if ( (_value + 1 ) >= (_geom->Size()[0]*(_geom->Size()[1] - 1)) - 1 )
			_valid = false;
	else if ( (_value + 2) % _geom->Size()[0] == 0 )
		_value += 3;
	else
		_value++;
	//std::cout<<_value<<" ";
}

BoundaryIterator::BoundaryIterator(const Geometry* geom) : Iterator(geom) {
	_boundary = 0;
	First();
}

void BoundaryIterator::SetBoundary(const index_t& boundary) {
	_boundary = boundary;
	First();
}

void BoundaryIterator::First() {
	_valid = true;

	if (_boundary == 0 )
	{
		_value = 0;
	}
	else if ( _boundary == 1 )
	{
		_value = 2*_geom->Size()[0] - 1;
	}
	else if ( _boundary == 2 )
	{
		_value = _geom->Size()[0]*_geom->Size()[1] - 1;
	}
	else if ( _boundary == 3 )
	{
		_value = _geom->Size()[0]*(_geom->Size()[1] - 2);
	}
}

void BoundaryIterator::Next() {
	if (_boundary == 0 )
	{
		_value++;
		if ( _value >= _geom->Size()[0] )
		{
			_valid = false;
		}
	}
	else if ( _boundary == 1 )
	{
		_value += _geom->Size()[0];
		if ( _value > _geom->Size()[0]*(_geom->Size()[1] - 1) )
		{
			_valid = false;
		}
	}
	else if ( _boundary == 2 )
	{
		_value--;
		if ( _value < _geom->Size()[0]*(_geom->Size()[1] - 1) )
		{
			_valid = false;
		}
	}
	else if ( _boundary == 3 )
	{
		_value -= _geom->Size()[0];
		if ( _value < _geom->Size()[0] )
		{
			_valid = false;
		}
	}
}
