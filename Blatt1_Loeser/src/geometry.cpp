#include "geometry.hpp"

#include <fstream>
#include <iostream>

using namespace std;

Geometry::Geometry() {
	_velocity[0] = 1;
	_velocity[1] = 1;
	_pressure = 1;

	_size[0] = 128;
	_size[1] = 128;

	_length[0] = 1;
	_length[1] = 1;

	_h[0] = _length[0]/_size[0];
	_h[1] = _length[1]/_size[1];
	cout << "Loaded default Variables:" << endl;
	PrintVariables();
}

void Geometry::Load(const char* file) {
	ifstream fileStream(file);

	if ( fileStream.fail() )
	{
		cout << "Couldn't load geometry from " << file << endl;
	}
	else
	{
		//Load values and compute mesh width
		fileStream >> _velocity[0];
		fileStream >> _velocity[1];
		fileStream >> _pressure;
		fileStream >> _size[0];
		fileStream >> _size[1];
		fileStream >> _length[0];
		fileStream >> _length[1];
		_h[0] = _length[0]/_size[0];
		_h[1] = _length[1]/_size[1];
		cout << "Loaded geometry from " << file << ":" << endl;
		PrintVariables();
	}


}
void Geometry::PrintVariables(){
	cout << "vel_x: " << _velocity[0] << endl;
	cout << "vel_y: " << _velocity[1] << endl;
	cout << "p: " << _pressure << endl;
	cout << "size_x: " << _size[0] << endl;
	cout << "size_y: " << _size[1] << endl;
	cout << "length_x: " << _length[0] << endl;
	cout << "length_y " << _length[1] << endl;
	cout << "h_x: " << _h[0] << endl;
	cout << "h_y: " << _h[1] << endl << endl;

}

const multi_index_t& Geometry::Size() const {
	return _size;
}

const multi_real_t& Geometry::Length() const {
	return _length;
}

const multi_real_t& Geometry::Mesh() const {
	return _h;
}

void Geometry::Update_U(Grid* u) const {
}

void Geometry::Update_V(Grid* v) const {
}

void Geometry::Update_P(Grid* p) const {
}
