#include "grid.hpp"

#include <iostream>
#include <cmath>

using namespace std;

Grid::Grid(const Geometry* geom) {
	_data = new real_t[geom->Size()[0]*geom->Size()[1]];
	_geom = geom;
	multi_real_t offset;
	offset[0] = 0;
	offset[1] = 0;
	_offset = offset;
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
	for (index_t i = 0; i< _geom->Size()[0]*_geom->Size()[1];i++){
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
	//TODO do proper interpolation
	index_t cell_x = (index_t)((pos[0] + _offset[0])/_geom->Mesh()[0]);
	index_t cell_y = (index_t)((pos[1] + _offset[1])/_geom->Mesh()[1]);

	return _data[cell_y*_geom->Size()[0] + cell_x];
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
	return (_data[it.Top()]-_data[it])/_geom->Mesh()[1];
}

real_t Grid::dxx(const Iterator& it) const {
	return (_data[it.Right()]-2*_data[it]+_data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
}

real_t Grid::dyy(const Iterator& it) const {
	return (_data[it.Top()]-2*_data[it]+_data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
}

real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const {

	real_t A = (_data[it]+_data[it.Right()])/2;
	real_t B = (_data[it]+_data[it.Left()])/2;
	return (( A*A - B*B)  +  alpha*(  abs(A)*A-abs(B)*B    )         )/_geom->Mesh()[0];
}

real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha,
		const Grid* v) const {
	real_t A = (v->Cell(it)+v->Cell(it.Right()))/2;
	real_t B = (_data[it]+_data[it.Top()])/2;
	//std::cout<<"B: "<<B;
	real_t C = (v->Cell(it.Down())+ v->Cell(it.Down().Right()))/2;
	real_t D = (_data[it]+_data[it.Down()] )/2;
	return ( A*B - C * D) + alpha*( abs(A)*B-abs(C)*D )  /_geom->Mesh()[1];
}

real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha,
		const Grid* u) const {
	real_t A = (u->Cell(it)+u->Cell(it.Top()))/2;
	real_t B = (_data[it]+_data[it.Right()])/2;
	real_t C = (u->Cell(it.Left())+ u->Cell(it.Left().Top()))/2;
	real_t D = (_data[it]+_data[it.Left()] )/2;
	return ( A*B - C * D) + alpha*( abs(A)*B-abs(C)*D )  /_geom->Mesh()[0];

}

real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
	real_t A = (_data[it]+_data[it.Top()])/2;
	real_t B = (_data[it]+_data[it.Down()])/2;
	return (( A*A - B*B)  +  alpha*(  abs(A)*A-abs(B)*B    )         )/_geom->Mesh()[1];
}

real_t Grid::Max() const {
	real_t max = _data[0];
	for (index_t i = 1; i< _geom->Size()[0]*_geom->Size()[1];i++){
		if (_data[i] >= max)
		max = _data[i];
	}
	return max;
}

real_t Grid::Min() const {
	real_t min = _data[0];
	for (index_t i = 1; i< _geom->Size()[0]*_geom->Size()[1];i++){
		if (_data[i] <= min)
		min = _data[i];
	}
	return min;
}

real_t Grid::AbsMax() const {
	real_t max = this->Max();
	real_t min = this->Min();
	if ((max+min) > 0)
		return max;
	else
		return min;

}

real_t* Grid::Data() {
	return _data;
}
