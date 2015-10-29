#include "compute.hpp"

Compute::Compute(const Geometry *geom, const Parameter *param)
{
	_geom = geom;
	_param = param;
}
  /// Deletes all grids
Compute::~Compute()
{
	delete _p;
	delete _u;
	delete _v;

	delete _F;
	delete _G;
	delete _rhs;

	delete _tmp;
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

// F = u + dt * A
// G = v + dt * B
void Compute::MomentumEqu(const real_t& dt) {
	InteriorIterator it = InteriorIterator(_geom);
	while ( it.Valid() )
	{
		const real_t u = _u->Cell(it);
		const real_t v = _v->Cell(it);

		const real_t RE_inv = 1.0/_param->Re();
		const real_t alpha = _param->Alpha();

		real_t A = RE_inv * ( _u->dxx(it) + _u->dyy(it) ) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
		real_t B = RE_inv * ( _v->dxx(it) + _v->dyy(it) ) - _v->DC_vdv_y(it, alpha) - _v->DC_udv_x(it, alpha, _u);

		_F->Cell(it) = u + dt*A;
		_G->Cell(it) = v + dt*B;

		it.Next();
	}
}

void Compute::RHS(const real_t& dt) {
	InteriorIterator it = InteriorIterator(_geom);
	while ( it.Valid() )
	{
		real_t dFx = (_F->Cell(it.Right()) - _F->Cell(it))/_geom->Mesh()[0];
		real_t dGy = (_G->Cell(it.Top()) - _G->Cell(it))/_geom->Mesh()[1];

		_rhs->Cell(it) = (dFx + dGy)/dt;
		it.Next();
	}
}
