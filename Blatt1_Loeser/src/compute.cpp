#include "compute.hpp"
#include <cmath>

Compute::Compute(const Geometry *geom, const Parameter *param)
{
	_geom = geom;
	_param = param;

	_solver = new SOR(_geom,_param->Omega());

	_t = param->Dt();
	//TODO was ist dtlimit?b
	_dtlimit = param->Tau();

	_epslimit = param->Eps();
	multi_real_t u_offset;
	multi_real_t v_offset;

	u_offset[0] = geom->Mesh()[0]/2;
	u_offset[1] = 0;
	v_offset[1] = geom->Mesh()[1]/2;
	v_offset[0] = 0;

	_u = new Grid(geom, u_offset);
	_v = new Grid(geom, v_offset);
	_p = new Grid(geom);
	_F = new Grid(geom);
	_G = new Grid(geom);
	_rhs = new Grid(geom);

	_u->Initialize(0);
	_v->Initialize(0);
	_p->Initialize(0);

	_tmp = new Grid(geom);
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
	delete _solver;
}

void Compute::TimeStep(bool printInfo) {

	//Compute dt
	// beta in (0,1).
	real_t beta = 0.9;
	real_t dt = beta*std::fmax(_geom->Mesh()[0],_geom->Mesh()[1])*std::fmax(_u->AbsMax(),_v->AbsMax());

	//Compute boudary Values
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);

	//Compute F, G
	MomentumEqu(dt);

	// Compute RHS
	RHS(dt);

	//Compute p
	real_t res = 10000000;
	index_t i = 0;
	while(res>_param->Eps()&& i<_param->IterMax()){
		res=_solver->Cycle(_p,_rhs);
		_geom->Update_P(_p);
		i++;
	}

	//Compute u,v
	NewVelocities(dt);

}

const real_t& Compute::GetTime() const {
	return _t;
}

const Grid* Compute::GetU() const {
	return _u;
}

const Grid* Compute::GetV() const {
	return _v;
}

const Grid* Compute::GetP() const {
	return _p;
}

const Grid* Compute::GetRHS() const {
	return _rhs;
}

const Grid* Compute::GetVelocity() {
	Iterator it = Iterator(_geom);

	while ( it.Valid() )
	{
		_tmp->Cell(it) = sqrt(( _u->Cell(it)*_u->Cell(it) + _v->Cell(it)*_v->Cell(it)));
		it.Next();
	}

	return _tmp;
}

const Grid* Compute::GetVorticity() {
}

const Grid* Compute::GetStream() {
}

void Compute::NewVelocities(const real_t& dt) {
	InteriorIterator it = InteriorIterator(_geom);
	while ( it.Valid() )
	{
		const real_t f = _F->Cell(it);
		const real_t g = _G->Cell(it);

		_u->Cell(it) = f - dt * ( _p->dx_l(it) );
		_v->Cell(it) = g - dt * ( _p->dy_l(it) );

		it.Next();
	}
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
