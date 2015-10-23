#include "compute.hpp"

Compute::Compute(const Geometry *geom, const Parameter *param)
{
	//TODO
}
  /// Deletes all grids
Compute::~Compute()
{
	//TODO
}

void Compute::TimeStep(bool printInfo) {
}

const real_t& Compute::GetTime() const {
}

const Grid* Compute::GetU() const {
}

const Grid* Compute::GetV() const {
}

const Grid* Compute::GetP() const {
}

const Grid* Compute::GetRHS() const {
}

const Grid* Compute::GetVelocity() {
}

const Grid* Compute::GetVorticity() {
}

const Grid* Compute::GetStream() {
}

void Compute::NewVelocities(const real_t& dt) {
}

void Compute::MomentumEqu(const real_t& dt) {
}

void Compute::RHS(const real_t& dt) {
}
