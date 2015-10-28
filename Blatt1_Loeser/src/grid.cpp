#include "grid.hpp"

#include <iostream>
//#include <cmath>
using namespace std;

Grid::Grid(const Geometry* geom) {
	_data = new real_t[geom->Size()[0]*geom->Size()[1]];
	_geom = geom;
}

Grid::Grid(const Geometry* geom, const multi_real_t& offset) {
	_data = new real_t[geom->Size()[0]*geom->Size()[1]];
	_offset = offset;
	_geom = geom;
}

Grid::~Grid() {
	delete[] _data;
}

void Grid::Initialize(const real_t& value) {
	for (int i = 0; i<= _geom->Size()[0]*_geom->Size()[1];i++){
		_data[i] = value;
	}
}

real_t& Grid::Cell(const Iterator& it) {
	return _data[it];
}

const real_t& Grid::Cell(const Iterator& it) const {
	return _data[it];
}

real_t Grid::Interpolate(const multi_real_t& pos) const {

}

real_t Grid::dx_l(const Iterator& it) const {
	return (_data[it]-_data[it.Left()])/_geom->Mesh()[0];
}

real_t Grid::dx_r(const Iterator& it) const {
	return (_data[it.Right()]-_data[it])/_geom->Mesh()[0];
}

real_t Grid::dy_l(const Iterator& it) const {
	return (_data[it]-_data[it.Down()])/_geom->Mesh()[1];
}

real_t Grid::dy_r(const Iterator& it) const {
	(_data[it.Top()]-_data[it])/_geom->Mesh()[1];
}

real_t Grid::dxx(const Iterator& it) const {
	(_data[it.Right()]-2*_data[it]+_data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
}

real_t Grid::dyy(const Iterator& it) const {
	(_data[it.Top()]-2*_data[it]+_data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
}

real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const {
}

real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha,
		const Grid* v) const {
}

real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha,
		const Grid* u) const {
}

real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
}

real_t Grid::Max() const {
}

real_t Grid::Min() const {
}

real_t Grid::AbsMax() const {
}

real_t* Grid::Data() {
}
