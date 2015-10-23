#include "grid.hpp"

Grid::Grid(const Geometry* geom) {
}

Grid::Grid(const Geometry* geom, const multi_real_t& offset) {
}

Grid::~Grid() {
}

void Grid::Initialize(const real_t& value) {
}

real_t& Grid::Cell(const Iterator& it) {
}

const real_t& Grid::Cell(const Iterator& it) const {
}

real_t Grid::Interpolate(const multi_real_t& pos) const {
}

real_t Grid::dx_l(const Iterator& it) const {
}

real_t Grid::dx_r(const Iterator& it) const {
}

real_t Grid::dy_l(const Iterator& it) const {
}

real_t Grid::dy_r(const Iterator& it) const {
}

real_t Grid::dxx(const Iterator& it) const {
}

real_t Grid::dyy(const Iterator& it) const {
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
