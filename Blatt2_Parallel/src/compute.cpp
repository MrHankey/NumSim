/*
 * Copyright (C) 2015  Raphael Leiteriz, Sebastian Reuschen, Hamzeh Kraus
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "compute.hpp"
#include <cmath>
#include <iostream>

#include "geometry.hpp"
#include "grid.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include "solver.hpp"
#include "communicator.hpp"

using namespace std;

/// Creates a compute instance with given geometry and parameter
//  @param geom  get geometry
//  @param param get parameter
Compute::Compute(const Geometry *geom, const Parameter *param, const Communicator *comm) {
	// Initialize
	_geom   = geom;
	_param  = param;
	_comm = comm;
	//_solver = new RedOrBlackSOR(_geom,_param->Omega());
	_solver = new SOR(_geom,_param->Omega());

	// Set time steps
	_t = 0.0;
	_dtlimit = param->Dt();

	// Define offset
	_epslimit = param->Eps();
	multi_real_t u_offset;
	multi_real_t v_offset;
	multi_real_t p_offset;
	multi_real_t vort_offset;
	multi_real_t stream_offset;

	// Set offset
	u_offset[0] = geom->Mesh()[0];
	u_offset[1] = geom->Mesh()[0]/2.0;
	v_offset[0] = geom->Mesh()[0]/2.0;
	v_offset[1] = geom->Mesh()[1];
	p_offset[0] = geom->Mesh()[0]/2.0;
	p_offset[1] = geom->Mesh()[1]/2.0;

	vort_offset[0] = geom->Mesh()[0];
	vort_offset[1] = geom->Mesh()[1];

	stream_offset[0] = geom->Mesh()[0];
	stream_offset[1] = geom->Mesh()[1];

	// Set grid for all variables
	_u   = new Grid(geom, u_offset);
	_v   = new Grid(geom, v_offset);
	_p   = new Grid(geom, p_offset);
	_F   = new Grid(geom);
	_G   = new Grid(geom);
	_temp1 = new Grid(geom);
	_temp2 = new Grid(geom);
	_rhs = new Grid(geom);
	_tmp = new Grid(geom);
	particelTracingGrid = new Grid(geom);
	StreakLinesGrid = new Grid(geom);
	particelTracingGrid->Initialize(0);

	//Start Position
	particelTracing.push_back({0,0.6});
	StreakLines.push_back({0,0.7});

	particelTracing_bot.push_back({0,0.2});
	StreakLines_bot.push_back({0,0.2});


	_tmpVorticity = new Grid(geom, vort_offset);
	_tmpStream    = new Grid(geom, stream_offset);

	// Set values
	_u->Initialize(0);
	_v->Initialize(0);
	_p->Initialize(0);

	_pL = 0.05; /// kann man auch noch aus Datei einlesen!
	_pR = 0.0;

	_tmpVorticity->Initialize(0);
	_tmpStream->Initialize(0);

}

/// Destructor: deletes all grids
Compute::~Compute() {
	delete _p;
	delete _u;
	delete _v;

	delete _F;
	delete _G;
	delete _rhs;

	delete _tmpVorticity;
	delete _tmpStream;

	delete _tmp;
	delete _solver;
}

/// Execute one time step of the fluid simulation (with or without debug info)
//  @param printInfo  print information about current solver state
void Compute::TimeStep(bool printInfo) {
	// Compute boundary Values
	// eigentlich erst nach dt, aber dann geht der erste Zeitschritt zu lang.
	//_geom->Update_U(_u);
	//_geom->Update_V(_v);
	//_geom->Update_P(_p);

	_geom->Update_All(_p,_u,_v,_pL,_pR);


	// Compute dt
	real_t dt = _param->Tau()*std::fmax(_geom->Mesh()[0],_geom->Mesh()[1])/std::fmax(_u->AbsMax(),_v->AbsMax());
	real_t dt2 = _param->Tau()*_param->Re()/2* (_geom->Mesh()[1]*_geom->Mesh()[1]*_geom->Mesh()[0]*_geom->Mesh()[0]);
	dt2 = dt2/(_geom->Mesh()[1]*_geom->Mesh()[1]+_geom->Mesh()[0]*_geom->Mesh()[0]);
	dt = std::min(dt2,std::min(dt,_param->Dt()));

	if(_t>0.8){
		multi_real_t partTrace = particelTracing.back();
		if(partTrace[0]<=_geom->Length()[0]){
			partTrace[0] = partTrace[0]+dt*_u->Interpolate(partTrace);
			partTrace[1] = partTrace[1]+dt*_v->Interpolate(partTrace);
			particelTracing.push_back(partTrace);
			//cout<<"x: "<<partTrace[0]<<" y: "<<partTrace[1]<<endl;

			//particelTracingGrid->Initialize(0);
			index_t cell_x = (index_t)ceil((partTrace[0]/_geom->Mesh()[0]));
			index_t cell_y = (index_t)ceil((partTrace[1]/_geom->Mesh()[1]));
			Iterator it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
			particelTracingGrid->Cell(it) = 1;
		}

		partTrace = particelTracing_bot.back();
		if(partTrace[0]<=_geom->Length()[0]){
			partTrace[0] = partTrace[0]+dt*_u->Interpolate(partTrace);
			partTrace[1] = partTrace[1]+dt*_v->Interpolate(partTrace);
			particelTracing_bot.push_back(partTrace);
			//cout<<"x: "<<partTrace[0]<<" y: "<<partTrace[1]<<endl;

			//particelTracingGrid->Initialize(0);
			index_t cell_x = (index_t)ceil((partTrace[0]/_geom->Mesh()[0]));
			index_t cell_y = (index_t)ceil((partTrace[1]/_geom->Mesh()[1]));
			Iterator it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
			particelTracingGrid->Cell(it) = 2;
		}

	}

	multi_real_t startPosition =  StreakLines.back();
	StreakLinesGrid->Initialize(0);
	for (multi_real_t &pos : StreakLines)
	{
		if(pos[0]<=_geom->Length()[0]){

			pos[0] = pos[0]+dt*_u->Interpolate(pos);
			pos[1] = pos[1]+dt*_v->Interpolate(pos);
			index_t cell_x = (index_t)ceil((pos[0]/_geom->Mesh()[0]));
			index_t cell_y = (index_t)ceil((pos[1]/_geom->Mesh()[1]));
			Iterator it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
			StreakLinesGrid->Cell(it) = 1;
			//cout<<"x: "<<pos[0]<<" y: "<<pos[1]<<endl;
		}
	}
	StreakLines.push_back(startPosition);
	index_t cell_x = (index_t)ceil((startPosition[0]/_geom->Mesh()[0]));
	index_t cell_y = (index_t)ceil((startPosition[1]/_geom->Mesh()[1]));
	Iterator it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
	StreakLinesGrid->Cell(it) = 1;

	startPosition =  StreakLines_bot.back();
	for (multi_real_t &pos : StreakLines_bot)
	{
		if(pos[0]<=_geom->Length()[0]){

			pos[0] = pos[0]+dt*_u->Interpolate(pos);
			pos[1] = pos[1]+dt*_v->Interpolate(pos);
			index_t cell_x = (index_t)ceil((pos[0]/_geom->Mesh()[0]));
			index_t cell_y = (index_t)ceil((pos[1]/_geom->Mesh()[1]));
			Iterator it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
			StreakLinesGrid->Cell(it) = 2;
			//cout<<"x: "<<pos[0]<<" y: "<<pos[1]<<endl;
		}
	}
	StreakLines_bot.push_back(startPosition);
	cell_x = (index_t)ceil((startPosition[0]/_geom->Mesh()[0]));
	cell_y = (index_t)ceil((startPosition[1]/_geom->Mesh()[1]));
	it = Iterator(_geom, cell_y*_geom->Size()[0] + cell_x);
	StreakLinesGrid->Cell(it) = 1;


	//biggest dt of all is chosen.
	//dt = _comm->gatherMin(dt);

	// Compute F, G
	MomentumEqu(dt);
	//_comm->copyBoundary(_F);
	//_comm->copyBoundary(_G);
	//std::cin.ignore();

	// Compute RHS
	RHS(dt);

	// Compute p
	real_t total_res = 10000000;
	real_t local_res = 0.0;
	index_t i = 0;

	// Update boundary
	while(total_res >_param->Eps() && i < _param->IterMax() ) {
		/*bool even = _comm->EvenOdd();
		if (even){
			local_res = _solver->RedCycle(_p,_rhs);
			_comm->copyBoundary(_p);
			local_res = _solver->BlackCycle(_p,_rhs);
			_comm->copyBoundary(_p);
		}
		else{
			local_res = _solver->BlackCycle(_p,_rhs);
			_comm->copyBoundary(_p);

			local_res = _solver->RedCycle(_p,_rhs);
			_comm->copyBoundary(_p);

		}*/
		local_res =_solver->Cycle(_p, _rhs);


		//_geom->Update_P(_p);
		_geom->Update_All(_p,_temp1,_temp2,_pL,_pR);
		total_res = local_res; //_comm->gatherMax(local_res)/_comm->getSize();//Max statt Sum da jetzt max gemacht wird
		i++;


	}

	// Compute u,v
	NewVelocities(dt);

	//cout << "new vels done" << endl;

	// Send u,v
	//_comm->copyBoundary(_u);
	//_comm->copyBoundary(_v);

	// Next timestep
	_t += dt;


	if (_t >= 49.0)
	{
		InteriorIterator it = InteriorIterator(_geom);
		index_t count = 0;
		real_t sum = 0;
		while (it.Valid())
		{

			multi_index_t it_pos = it.Pos();
			multi_real_t sample_pos = {(it_pos[0]/(real_t)_geom->Size()[0])*_geom->Length()[0], (it_pos[1]/(real_t)_geom->Size()[1])*_geom->Length()[1] };

			real_t exact = -0.5*_param->Re()*(_pL-_pR)/_geom->Length()[0]*sample_pos[1]*(sample_pos[1] - _geom->Length()[1]);


			sum += fabs(_u->Interpolate(sample_pos) - exact);

			count++;
			it.Next();
		}

		sum /= count;
		cout << "L2 error: " << sum << endl;
	}

	// Print info
	if (printInfo) {
		cout << "t: " << _t << " dt: " << dt << " iter: " << i << "  \tres: " << std::scientific << total_res << "\t progress: " << std::fixed << _t/_param->Tend()*100 << "%" << endl;
		//cout<<"Theoretisches Maximum: "<<_param->Re()*(_pL-_pR)/_geom->Length()[0]*_geom->Length()[1]*_geom->Length()[1]/8.0<<" Unser maximum: "<<_u->AbsMax()<<endl;

		/*multi_real_t sample_pos_1 = {_geom->Length()[0]/2, _geom->Length()[1]/6 };
		multi_real_t sample_pos_2 = {_geom->Length()[0]/2, _geom->Length()[1]/3 };
		multi_real_t sample_pos_3 = {_geom->Length()[0]/2, _geom->Length()[1]/2 };

		real_t exact1 = -0.5*_param->Re()*(_pL-_pR)/_geom->Length()[0]*sample_pos_1[1]*(sample_pos_1[1] - _geom->Length()[1]);
		real_t exact2 = -0.5*_param->Re()*(_pL-_pR)/_geom->Length()[0]*sample_pos_2[1]*(sample_pos_2[1] - _geom->Length()[1]);
		real_t exact3 = -0.5*_param->Re()*(_pL-_pR)/_geom->Length()[0]*sample_pos_3[1]*(sample_pos_3[1] - _geom->Length()[1]);

		real_t diff1 = fabs(_u->Interpolate(sample_pos_1) - exact1);
		real_t diff2 = fabs(_u->Interpolate(sample_pos_2) - exact2);
		real_t diff3 = fabs(_u->Interpolate(sample_pos_3) - exact3);

		cout << "diffs_" << _geom->Size()[1] - 2 << " = [" << diff1 << ", " << diff2 << " ," << diff3 << "]" << endl;*/
	}

	//_u->Cell(it) = -10000;
}

/* Getter functions */
/// Returns the simulated time in total
const real_t& Compute::GetTime() const {return _t;}
/// Returns the pointer to U
const Grid*   Compute::GetU()    const {return _u;}
/// Returns the pointer to V
const Grid*   Compute::GetV()    const {return _v;}
/// Returns the pointer to P
const Grid*   Compute::GetP()    const {return _p;}
/// Returns the pointer to RHS
const Grid*   Compute::GetRHS()  const {return _rhs;}


/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity() {
	// Initialize
	Iterator it = Iterator(_geom);

	// Cycle through all cells
	while ( it.Valid() ) {
		_tmp->Cell(it) = sqrt(( _u->Cell(it)*_u->Cell(it) + _v->Cell(it)*_v->Cell(it)));
		it.Next();
	}

	return _tmp;
}

/// Computes and returns the vorticity
const Grid* Compute::GetVorticity() {
	// Initialize
	Iterator it = InteriorIterator(_geom);

	// Cycle through all cells
	if(false) {
		while(it.Valid()) {
			_tmpVorticity->Cell(it) = _u->dy_r(it) - _v->dx_r(it);
			it.Next();
		}
	} else {
		while(it.Valid()) {
			if     (_u->dy_r(it) - _v->dx_r(it)>0) _tmpVorticity->Cell(it) =  1;
			else if(_u->dy_r(it) - _v->dx_r(it)<0) _tmpVorticity->Cell(it) = -1;
			else								   _tmpVorticity->Cell(it) =  0;
			it.Next();
		}
		BoundaryIterator itb = BoundaryIterator(_geom);
		itb.SetBoundary(0);

		while(itb.Valid()) {
			if(_tmpVorticity->Cell(itb.Top()) == 1 /*&& _tmpVorticity->Cell(itb.Top().Left()) == -1*/) {
				_tmpVorticity->Cell(itb.Top()) = 2;
				std::cout<<"d_breakoff: "<<(itb.Value()*_geom->Mesh()[0])<<endl;
				break;
			}
			itb.Next();
		}
	}
	return _tmpVorticity;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream() {
	// Initialize
	Iterator it = InteriorIterator(_geom);

	multi_index_t screenSize;
	screenSize[0] = (_geom->TotalSize()[0]-2)/(_geom->Size()[0]-2);
	screenSize[1] = (_geom->TotalSize()[1]-2)/(_geom->Size()[1]-2);

	for(int i = 0;i<screenSize[0]+screenSize[1];i++){

		// Set first value to zero
		it.First();
		_tmpStream->Cell(it.Left()) = 0;
		//_comm->copyBoundary(_tmpStream);
		//it.Next();

		// Cycle through all cells
		while (it.Valid()){
			if (it.Pos()[1] == 0){
				_tmpStream->Cell(it) = _tmpStream->Cell(it.Left()) - _v->Cell(it)*_geom->Mesh()[0];
			} else {
				_tmpStream->Cell(it) = _tmpStream->Cell(it.Down()) + _u->Cell(it)*_geom->Mesh()[1];
			}
			it.Next();
		}
	}
	// TODO communicate boundary

	return _tmpStream;
}

/// Compute the new velocites u,v
//  @param dt  get timestep
void Compute::NewVelocities(const real_t& dt) {
	// Initialize interior iterator
	InteriorIterator it = InteriorIterator(_geom);

	// Cycle through all interior cells
	while ( it.Valid() ) {
		// Get temporary velocities
		const real_t f = _F->Cell(it);
		const real_t g = _G->Cell(it);

		// Calculate new velocities
		_u->Cell(it) = f - dt * ( _p->dx_r(it) );
		_v->Cell(it) = g - dt * ( _p->dy_r(it) );

		// Next cell
		it.Next();
	}
}


/* F = u + dt * A */
/* G = v + dt * B */
/// Compute the temporary velocities F,G
//  @param dt  get timestep
void Compute::MomentumEqu(const real_t& dt) {
	// Initialize interior iterator
	Iterator it = Iterator(_geom);
	it.First();

	//updaten oder nicht updaten?
	//_geom->Update_All(_temp1,_F,_G,_pL,_pR);

	// Cycle through all interior cells
	while ( it.Valid() ) {
		const real_t u = _u->Cell(it);
		const real_t v = _v->Cell(it);

		// Update new temporary velocities
		if ( _geom->_b->Cell(it) == 0)
		{
			// Get current velocities
			//const real_t u = _u->Cell(it);
			//const real_t v = _v->Cell(it);

			// Get parameter
			const real_t RE_inv = 1.0/_param->Re();
			const real_t alpha = _param->Alpha();

			// Calculate temporary velocietes
			real_t A = RE_inv * ( _u->dxx(it) + _u->dyy(it) ) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
			real_t B = RE_inv * ( _v->dxx(it) + _v->dyy(it) ) - _v->DC_vdv_y(it, alpha) - _v->DC_udv_x(it, alpha, _u);


			_F->Cell(it) = u + dt*A;
			_G->Cell(it) = v + dt*B;
		}
		else if (_geom->_b->Cell(it) == 1 )
		{
			_G->Cell(it) = v;
		}
		else //( _geom->_b->Cell(it) != 0 )
		{
			_F->Cell(it) = u;
		}
		//cout<<_F->Cell(it)<<" das war bei it = "<<it.Value()<<endl;

		// Next cell
		it.Next();
	}
	//_geom->Update_All(_temp1,_F,_G,_pL,_pR);
	//it.First();
	//cout<<_F->Cell(it)<<" das war bei it after: = "<<it<<endl;
	//_geom->Update_U(_F);
	//_geom->Update_V(_G);
}

/// Compute the RHS of the poisson equation
//  @param dt  get timestep
void Compute::RHS(const real_t& dt) {
	// Initialize interior iterator
	InteriorIterator it = InteriorIterator(_geom);

	// Cycle through all interior cells
	while ( it.Valid() ) {
		//real_t dFx = (_F->Cell(it.Right()) - _F->Cell(it))/_geom->Mesh()[0];
		//real_t dGy = (_G->Cell(it.Top())   - _G->Cell(it))/_geom->Mesh()[1];

		// Calculate RHS for both temporary velocities
		real_t dFx = (_F->Cell(it) - _F->Cell(it.Left()))/_geom->Mesh()[0];
		real_t dGy = (_G->Cell(it) - _G->Cell(it.Down()))/_geom->Mesh()[1];

		// Update new RHS
		_rhs->Cell(it) = (dFx + dGy)/dt;
		//std::cout<<_rhs->Cell(it)<<" "<<endl;

		// Next cell
		it.Next();
	}
}
