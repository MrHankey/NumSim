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
using namespace cl;

/// Creates a compute instance with given geometry and parameter
//  @param geom  get geometry
//  @param param get parameter
Compute::Compute(const Geometry *geom, const Parameter *param, OCLManager* oclmanager) {
	// Initialize
	_geom   = geom;
	_param  = param;
	_oclmanager = oclmanager;
	//_solver = new RedOrBlackSOR(_geom,_param->Omega());
	//_solver = new SOR(_geom,_param->Omega());
	//_solver = new JacobiOCL(_geom);
	_solver = new SOROCL(_geom, _param->Omega(), _oclmanager);

	_solver_time = 0.0;
	_time_buf = 0.0f;
	_time_buf_read = 0.0f;
	_time_kernel = 0.0f;
	_time_res = 0.0f;

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
	_rhs = new Grid(geom);
	_tmp = new Grid(geom);
	_T = new Grid(geom, p_offset);
	_locRes = new Grid(geom, p_offset);

	_tmpVorticity = new Grid(geom, vort_offset);
	_tmpStream    = new Grid(geom, stream_offset);

	// Set values
	_u->Initialize(0);
	_v->Initialize(0);
	_p->Initialize(0);
	_F->Initialize(0);
	_G->Initialize(0);
	_rhs->Initialize(0);
	_T->Initialize(0);
	_locRes->Initialize(0);

	_tmpVorticity->Initialize(0);
	_tmpStream->Initialize(0);

}

/// Destructor: deletes all grids
Compute::~Compute() {
	delete _p;
	delete _u;
	delete _v;
	delete _T;

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

	static real_t timing_dt = 0.0f;
	clock_t begin;
	begin = clock();

	// Compute dt
	real_t h_x = _geom->Mesh()[0];
	real_t h_y = _geom->Mesh()[1];
	real_t hs_x = h_x*h_x;
	real_t hs_y = h_y*h_y;
	// CFL
	real_t dt = fmax(h_x,h_y)/_oclmanager->ReduceMaxVelocity();
	// Mesh width and RE
	real_t dt2 = _param->Re()/2* (hs_y*hs_x);
	dt2 = dt2/(hs_y+hs_x);
	dt = std::min(dt2,std::min(dt,_param->Dt()));
#ifdef HEAT_COUPLING
	// Temp restriction
	real_t dt3 = (_param->Re()*_param->Pr()*0.5f)*(1.0f/(1.0f/hs_x + 1.0f/hs_y));
	//real_t dt3 = 0.5 * _param->Re() * _param->Pr() / ( 1/hs_x + 1/hs_y );
	dt = std::min(dt, dt3);
#endif

	dt *= _param->Tau();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//_solver_time = _solver->_gpu_time;
	timing_dt += elapsed_secs;

#ifdef HEAT_COUPLING
	Temperature(dt);
#endif

	static real_t timing_meq = 0.0f;
	begin = clock();
	// Compute F, G
	MomentumEqu(dt);
	_oclmanager->_queue.finish();
	//_comm->copyBoundary(_F);
	//_comm->copyBoundary(_G);
	//std::cin.ignore();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	timing_meq += elapsed_secs;

	static real_t timing_rhs = 0.0f;
	begin = clock();
	// Compute RHS
	RHS(dt);
	_oclmanager->_queue.finish();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	timing_rhs += elapsed_secs;

	// Compute p
	real_t total_res = 10000000;
	index_t i = 0;
	index_t numBlockCycles = 10;

	begin = clock();

	while(total_res >_param->Eps() && i < _param->IterMax() ) {
		total_res = _solver->Cycle(numBlockCycles);
		i += numBlockCycles;
	}

	_oclmanager->_queue.finish();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	_solver_time += elapsed_secs;


	_time_buf = _solver->_time_buffer;
	_time_buf_read = _solver->_time_buffer_read;
	_time_kernel = _solver->_time_kernel;
	_time_res = _solver->_time_res;

	static real_t timing_nv = 0.0f;
	begin = clock();
	// Compute u,v
	NewVelocities(dt);
	_oclmanager->_queue.finish();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	timing_nv += elapsed_secs;

	//cout << "new vels done" << endl;

	// Send u,v
	//_comm->copyBoundary(_u);
	//_comm->copyBoundary(_v);

	// Next timestep
	_t += dt;

	static real_t read_time = 0.0f;

	begin = clock();
	// Print info
	if (printInfo) {
		//_oclmanager->_queue.finish();
		cout << "t: " << _t << " dt: " << dt << "  \tres: " << std::scientific << total_res << "\t progress: " << std::fixed << _t/_param->Tend()*100 << "%" << " IterCount: " << i << endl;
		if ( _t >= _param->Tend())
		{
			cout << "Timestep timings:" << endl;
			cout << "dt_time: \t" << timing_dt << " \nmeq: \t\t" << timing_meq << " \nrhs: \t\t" << timing_rhs << " \nsolver: \t" << _solver_time << " \nnv: \t\t" << timing_nv << " \nread: \t\t" << read_time << endl << endl;
		}
	}
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	read_time += elapsed_secs;

}

/* Getter functions */
/// Returns the simulated time in total
const real_t& Compute::GetTime() const {return _t;}
/// Returns the pointer to U
const Grid*   Compute::GetU()    const {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_u, CL_TRUE, 0, gridSize * sizeof(real_t), _u->_data);
	return _u;
}
/// Returns the pointer to V
const Grid*   Compute::GetV()    const {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_v, CL_TRUE, 0, gridSize * sizeof(real_t), _v->_data);
	return _v;
}
/// Returns the pointer to P
const Grid*   Compute::GetP()    const {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_p, CL_TRUE, 0, gridSize * sizeof(real_t), _p->_data);
	return _p;
}
/// Returns the pointer to RHS
const Grid*   Compute::GetRHS()  const {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_rhs, CL_TRUE, 0, gridSize * sizeof(real_t), _rhs->_data);
	return _rhs;
}
/// Returns the pointer to T
const Grid*   Compute::GetT()  const {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_T, CL_TRUE, 0, gridSize * sizeof(real_t), _T->_data);
	return _T;
}

const Grid*	  Compute::GetRes() const
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_locRes, CL_TRUE, 0, gridSize * sizeof(real_t), _locRes->_data);
	return _locRes;
}


/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity() {

	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_u, CL_TRUE, 0, gridSize * sizeof(real_t), _u->_data);
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_v, CL_TRUE, 0, gridSize * sizeof(real_t), _v->_data);
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
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_u, CL_TRUE, 0, gridSize * sizeof(real_t), _u->_data);
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_v, CL_TRUE, 0, gridSize * sizeof(real_t), _v->_data);
	// Initialize
	Iterator it = InteriorIterator(_geom);

	// Cycle through all cells
	while(it.Valid()) {
		_tmpVorticity->Cell(it) = _u->dy_r(it) - _v->dx_r(it);
		it.Next();
	}

	return _tmpVorticity;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream() {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_u, CL_TRUE, 0, gridSize * sizeof(real_t), _u->_data);
	_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_v, CL_TRUE, 0, gridSize * sizeof(real_t), _v->_data);
	// Initialize
	Iterator it = InteriorIterator(_geom);

	multi_index_t screenSize;
	screenSize[0] = (_geom->Size()[0]-2)/(_geom->Size()[0]-2);
	screenSize[1] = (_geom->Size()[1]-2)/(_geom->Size()[1]-2);

	for(index_t i = 0;i<screenSize[0]+screenSize[1];i++){

		// Set first value to zero
		it.First();
		_tmpStream->Cell(it.Left()) = 0;

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

	index_t localSize = 16;
	real_t dt_var = dt;
	real_t fVelocity = _geom->Vel()[0];


	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_F, CL_TRUE, 0, gridSize * sizeof(real_t), _F->_data);
	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_G, CL_TRUE, 0, gridSize * sizeof(real_t), _G->_data);

	Buffer clDT = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &dt_var);
	Buffer clVelocity = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fVelocity);
	// Set arguments to kernel
	checkErr(_oclmanager->_kernel_newvel.setArg(0, _oclmanager->_F), "setArg0");
	checkErr(_oclmanager->_kernel_newvel.setArg(1, _oclmanager->_G), "setArg1");
	checkErr(_oclmanager->_kernel_newvel.setArg(2, _oclmanager->_p), "setArg2");
	checkErr(_oclmanager->_kernel_newvel.setArg(3, _oclmanager->_u), "setArg3");
	checkErr(_oclmanager->_kernel_newvel.setArg(4, _oclmanager->_v), "setArg4");
	checkErr(_oclmanager->_kernel_newvel.setArg(5, _oclmanager->_h_inv), "setArg5");
	checkErr(_oclmanager->_kernel_newvel.setArg(6, clDT), "setArg6");
	checkErr(_oclmanager->_kernel_newvel.setArg(7, clVelocity), "setArg7");

	NDRange global((_geom->Size()[0] - 2), (_geom->Size()[1] - 2));
	NDRange local(localSize,localSize);
	checkErr(_oclmanager->_queue.enqueueNDRangeKernel(_oclmanager->_kernel_newvel, NullRange, global, local), "enqueueNDRangeKernel");

	// Initialize interior iterator
	/*InteriorIterator it = InteriorIterator(_geom);

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
	}*/
}


/* F = u + dt * A */
/* G = v + dt * B */
/// Compute the temporary velocities F,G
//  @param dt  get timestep
void Compute::MomentumEqu(const real_t& dt) {

	index_t localSize = 16;
	real_t dt_var = dt;
	real_t RE_inv = 1.0f/_param->Re();
	real_t fGravity = _param->Gravity();
	real_t fBeta = _param->Beta();
	real_t fVelocity = _geom->Vel()[0];

	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_u, CL_TRUE, 0, gridSize * sizeof(real_t), _u->_data);
	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_v, CL_TRUE, 0, gridSize * sizeof(real_t), _v->_data);

	Buffer clDT = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &dt_var);
	Buffer clRE_inv = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &RE_inv);
	Buffer clGravity = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fGravity);
	Buffer clBeta = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fBeta);
	Buffer clVelocity = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fVelocity);
	// Set arguments to kernel
	checkErr(_oclmanager->_kernel_momentumeq.setArg(0, _oclmanager->_F), "setArg0 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(1, _oclmanager->_G), "setArg1 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(2, _oclmanager->_u), "setArg2 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(3, _oclmanager->_v), "setArg3 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(4, _oclmanager->_h_inv), "setArg4 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(5, _oclmanager->_hs_inv), "setArg5 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(6, clDT), "setArg6 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(7, clRE_inv), "setArg7 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(8, _oclmanager->_T), "setArg8 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(9, clGravity), "setArg9 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(10, clBeta), "setArg10 meq");
	checkErr(_oclmanager->_kernel_momentumeq.setArg(11, clVelocity), "setArg11 meq");

	NDRange global((_geom->Size()[0] - 2), (_geom->Size()[1] - 2));
	NDRange local(localSize,localSize);
	checkErr(_oclmanager->_queue.enqueueNDRangeKernel(_oclmanager->_kernel_momentumeq, NullRange, global, local), "enqueueNDRangeKernel meq");

	//_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_F, CL_TRUE, 0, gridSize * sizeof(real_t), _F->_data);
	//_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_G, CL_TRUE, 0, gridSize * sizeof(real_t), _G->_data);
	//_oclmanager->_queue.finish();

	// Initialize interior iterator
	//InteriorIterator it = InteriorIterator(_geom);

	// Cycle through all interior cells
	/*while ( it.Valid() ) {
		// Get current velocities
		const real_t u = _u->Cell(it);
		const real_t v = _v->Cell(it);

		// Get parameter
		const real_t RE_inv = 1.0/_param->Re();
		const real_t alpha = _param->Alpha();

		// Calculate temporary velocietes
		real_t A = RE_inv * ( _u->dxx(it) + _u->dyy(it) ) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
		real_t B = RE_inv * ( _v->dxx(it) + _v->dyy(it) ) - _v->DC_vdv_y(it, alpha) - _v->DC_udv_x(it, alpha, _u);

		// Update new temporary velocities
		_F->Cell(it) = u + dt*A;
		_G->Cell(it) = v + dt*B;

		// Next cell
		it.Next();
	}*/

	//_geom->Update_U(_F);
	//_geom->Update_V(_G);
}

/// Compute the RHS of the poisson equation
//  @param dt  get timestep
void Compute::RHS(const real_t& dt) {

	index_t localSize = 16;
	real_t dt_inv = 1.0f/dt;


	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_F, CL_TRUE, 0, gridSize * sizeof(real_t), _F->_data);
	//_oclmanager->_queue.enqueueWriteBuffer(_oclmanager->_G, CL_TRUE, 0, gridSize * sizeof(real_t), _G->_data);

	Buffer clDTInv = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &dt_inv);
	// Set arguments to kernel
	checkErr(_oclmanager->_kernel_rhs.setArg(0, _oclmanager->_F), "setArg0 rhs");
	checkErr(_oclmanager->_kernel_rhs.setArg(1, _oclmanager->_G), "setArg1 rhs");
	checkErr(_oclmanager->_kernel_rhs.setArg(2, _oclmanager->_rhs), "setArg2 rhs");
	checkErr(_oclmanager->_kernel_rhs.setArg(3, _oclmanager->_h_inv), "setArg3 rhs");
	checkErr(_oclmanager->_kernel_rhs.setArg(4, clDTInv), "setArg4 rhs");

	NDRange global((_geom->Size()[0] - 2), (_geom->Size()[1] - 2));
	NDRange local(localSize,localSize);
	checkErr(_oclmanager->_queue.enqueueNDRangeKernel(_oclmanager->_kernel_rhs, NullRange, global, local), "enqueueNDRangeKernel");

	//_oclmanager->_queue.enqueueReadBuffer(_oclmanager->_rhs, CL_TRUE, 0, gridSize * sizeof(real_t), _rhs->_data);
	//_oclmanager->_queue.finish();

	// Initialize interior iterator
	/*InteriorIterator it = InteriorIterator(_geom);

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
	}*/
}

/// Compute the Temperature
//  @param dt  get timestep
void Compute::Temperature(const real_t& dt) {

	index_t localSize = 16;
	real_t dt_var = dt;

	real_t RE_inv = 1.0f/_param->Re();
	real_t PR_inv = 1.0f/_param->Pr();
	real_t fTemp = _param->T();
	real_t fQ = _param->Q();
	real_t fGamma = _param->Gamma();

	Buffer clDT = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &dt_var);
	Buffer clRE_inv = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &RE_inv);
	Buffer clPR_inv = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &PR_inv);
	Buffer clTemp = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fTemp);
	Buffer clQ = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fQ);
	Buffer clGamma = Buffer(_oclmanager->_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &fGamma);
	// Set arguments to kernel
	checkErr(_oclmanager->_kernel_T.setArg(0, _oclmanager->_T), "setArg0 T");
	checkErr(_oclmanager->_kernel_T.setArg(1, _oclmanager->_u), "setArg1 T");
	checkErr(_oclmanager->_kernel_T.setArg(2, _oclmanager->_v), "setArg2 T");
	checkErr(_oclmanager->_kernel_T.setArg(3, _oclmanager->_h_inv), "setArg3 temp");
	checkErr(_oclmanager->_kernel_T.setArg(4, _oclmanager->_hs_inv), "setArg4 temp");
	checkErr(_oclmanager->_kernel_T.setArg(5, clDT), "setArg5 temp");
	checkErr(_oclmanager->_kernel_T.setArg(6, clRE_inv), "setArg6 temp");
	checkErr(_oclmanager->_kernel_T.setArg(7, clPR_inv), "setArg7 temp");
	checkErr(_oclmanager->_kernel_T.setArg(8, clTemp), "setArg8 temp");
	checkErr(_oclmanager->_kernel_T.setArg(9, clQ), "setArg9 temp");
	checkErr(_oclmanager->_kernel_T.setArg(10, clGamma), "setArg10 temp");

	NDRange global((_geom->Size()[0] - 2), (_geom->Size()[1] - 2));
	NDRange local(localSize,localSize);
	checkErr(_oclmanager->_queue.enqueueNDRangeKernel(_oclmanager->_kernel_T, NullRange, global, local), "enqueueNDRangeKernel temp");


}
