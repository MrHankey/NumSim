#include "solver.hpp"

Solver::Solver(const Geometry* geom) {
}

Solver::~Solver() {
}

real_t Solver::localRes(const Iterator& it, const Grid* grid,
		const Grid* rhs) const {
}

SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
}

SOR::~SOR() {
}

real_t SOR::Cycle(Grid* grid, const Grid* rhs) const {
}
