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

#include "solver.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

#include <cmath>
#include <utility>
#include <fstream>

using namespace cl;
using namespace std;

/// Constructor of the abstract Solver class
//  @param geom  get field geometry
Solver::Solver(const Geometry* geom) {
	_geom = geom;
}

/// Destructor of the Solver Class
Solver::~Solver() {}

/// Returns the residual at [it] for the pressure-Poisson equation
//  @param  it        calls iterator
//  @param  grid      get grid values
//  @return localRes  return residual
real_t Solver::localRes(const Iterator& it, const Grid* grid,
		const Grid* rhs) const {
	return fabs(grid->dxx(it) + grid->dyy(it) - rhs->Cell(it));
}

/* SOR solver */
/// Concrete SOR solver
//  @param geom   get geometry
//  @param omega  get scaling factor
SOR::SOR(const Geometry* geom, const real_t& omega) : Solver(geom) {
	_geom  = geom;
	_omega = omega;
}

/// Destructor
SOR::~SOR() {}


/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t SOR::Cycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy, norm;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	// Cycle through all cells
	while(it.Valid()) {
		// Initialize
		n++;
		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		// define variables
		pij  	= grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t 	= grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r 	= grid->Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		// Save new p_ij
		grid->Cell(it)  = (1-_omega)*pij + _omega * norm * corr;
		// Calculate and summ residual
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		// Next cell
		it.Next();
	}

	// Norm residual
	return res/n;
}

//---------------------------------------------------------------------------------------------------

/* SOR solver */
/// Concrete SOR solver
//  @param geom   get geometry
//  @param omega  get scaling factor
JacobiOCL::JacobiOCL( const Geometry* geom) {

	_geom  = geom;

	vector<Platform> platforms;
	Platform::get(&platforms);

	// Select the default platform and create a context using this platform and the GPU
	cl_context_properties cps[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)(platforms[0])(),
		0
	};
	_context = Context( CL_DEVICE_TYPE_ALL, cps);

	// Get a list of devices on this platform
	//vector<Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	platforms[0].getDevices(CL_DEVICE_TYPE_ALL, &_all_devices);
	if(_all_devices.size()==0){
		std::cout<<" No devices found. Check OpenCL installation!\n";
		exit(1);
	}

	// Create a command queue and use the first device
	_queue = CommandQueue(_context, _all_devices[0]);

	//cout << "Using device: " << _all_devices[0].getInfo<CL_DEVICE_NAME>() << endl;

	// Read source file
	std::ifstream sourceFile("poisson.cl");
	std::string sourceCode(
		std::istreambuf_iterator<char>(sourceFile),
		(std::istreambuf_iterator<char>()));
	Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	// Make program of the source code in the context
	_program = Program(_context, source);

	// Build program for these specific devices
	_program.build(_all_devices);
	std::cout<<" Error building: "<<_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_all_devices[0])<<"\n";
}

/// Destructor
JacobiOCL::~JacobiOCL() {}


/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t JacobiOCL::Cycle(Grid* grid, const Grid* rhs) {

	real_t dx   = _geom->Mesh()[0];
	real_t n = _geom->Size()[0]*_geom->Size()[1];

	/*real_t *A = new real_t[(index_t)n];
	real_t *B = new real_t[(index_t)n];
	for(int i = 0; i < n; i++) {
		A[i] = (double)i;
		B[i] = (double)(n - i);
	}*/

	// Make kernel
	_kernel = Kernel(_program, "poisson_jacobi");

	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	// Create memory buffers
	_bufOld = Buffer(_context, CL_MEM_READ_ONLY, gridSize * sizeof(real_t));
	_bufRHS = Buffer(_context, CL_MEM_READ_ONLY, gridSize * sizeof(real_t));
	_bufNew = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t));
	Buffer clDx = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &dx);

#ifdef __CL_ENABLE_EXCEPTIONS
	try	{
#endif

		// Copy lists A and B to the memory buffers
		_queue.enqueueWriteBuffer(_bufOld, CL_TRUE, 0, (index_t)n * sizeof(real_t), grid->_data);
		_queue.enqueueWriteBuffer(_bufRHS, CL_TRUE, 0, (index_t)n * sizeof(real_t), rhs->_data);
		_queue.enqueueWriteBuffer(_bufNew, CL_TRUE, 0, (index_t)n * sizeof(real_t), grid->_data);

		// Set arguments to kernel
		_kernel.setArg(0, _bufOld);
		_kernel.setArg(1, _bufRHS);
		_kernel.setArg(2, _bufNew);
		//_kernel.setArg(3, sizeof(clDx), &clDx);
		_kernel.setArg(3, clDx);

		// Run the kernel on specific ND range
		NDRange global(_geom->Size()[0] - 2, _geom->Size()[1] - 2);
		NDRange local(1,1);
		_queue.enqueueNDRangeKernel(_kernel, NullRange, global, local);
		//_queue.enqueue


		/*Grid newGrid = Grid(_geom);
		newGrid.Initialize(0.0);*/
		_queue.enqueueReadBuffer(_bufNew, CL_TRUE, 0, n * sizeof(real_t), grid->_data);
#ifdef __CL_ENABLE_EXCEPTIONS
	} catch(Error error) {
	   std::cout << "Error initializing OpenCL: " << error.what() << "(" << error.err() << ")" << std::endl;
	   exit(1);
	}
#endif

	_queue.finish();
	InteriorIterator it = InteriorIterator(_geom);
	int num = 0;
	real_t res = 0;
	while(it.Valid()) {
		// Calculate and summ residual
		real_t dxx = grid->dxx(it);
		real_t dyy = grid->dyy(it);
		real_t fRhs = rhs->Cell(it);
		res += fabs( dxx + dyy - fRhs);
		it.Next();
		num++;
	}

	// Norm residual
	return res/num;

	/*for(int i = 0; i < LIST_SIZE; i ++)
		 std::cout << A[i] << " + " << B[i] << " = " << C[i] << std::endl; */
}

//---------------------------------------------------------------------------------------------------

/* SOR solver */
/// Concrete SOR solver
//  @param geom   get geometry
//  @param omega  get scaling factor
RedOrBlackSOR::RedOrBlackSOR(const Geometry* geom, const real_t& omega) : SOR(geom,omega) {
	_geom  = geom;
	_omega = omega;
}

/// Destructor
RedOrBlackSOR::~RedOrBlackSOR() {}

/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t RedOrBlackSOR::RedCycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy, norm;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	// Cycle through all cells
	while(it.Valid()) {
		// Initialize
		n++;
		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		// define variables
		pij  	= grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t 	= grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r 	= grid->Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		// Save new p_ij
		grid->Cell(it)  = (1-_omega)*pij + _omega * norm * corr;
		// Calculate and summ residual
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		// Next cell
		it.Next();
		it.Next();
	}

	// Norm residual
	return res/n;
}


/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t RedOrBlackSOR::BlackCycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy, norm;
	it.Next();
	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	// Cycle through all cells
	while(it.Valid()) {
		// Initialize
		n++;
		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		// define variables
		pij  	= grid->Cell(it);
		pij_d   = grid->Cell(it.Down());
		pij_t 	= grid->Cell(it.Top());
		pij_l   = grid->Cell(it.Left());
		pij_r 	= grid->Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		// Save new p_ij
		grid->Cell(it)  = (1-_omega)*pij + _omega * norm * corr;
		// Calculate and summ residual
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		// Next cell
		it.Next();
		it.Next();
	}

	// Norm residual
	return res/n;
}

Jacobi::Jacobi(const Geometry* geom) : Solver(geom) {
	_geom = geom;
}

Jacobi::~Jacobi() {
}

real_t Jacobi::Cycle(Grid* grid, const Grid* rhs) const {
	// Initialize
	InteriorIterator it = InteriorIterator(_geom);
	real_t n, res, dx, dy;

	n    = 0;
	res  = 0;
	dx   = _geom->Mesh()[0];
	dy   = _geom->Mesh()[1];
	//norm = 0.5*( (dx*dx*dy*dy)/(dx*dx+dy*dy) );

	Grid oldGrid = Grid(_geom);
	while(it.Valid()) {
		oldGrid.Cell(it) = grid->Cell(it);
		it.Next();

	}

    it.First();
	// Cycle through all cells
	while(it.Valid()) {
		// Initialize
		n++;
		real_t pij, pij_d, pij_l, pij_t, pij_r, A, B, corr;

		// define variables
		pij  	= oldGrid.Cell(it);
		pij_d   = oldGrid.Cell(it.Down());
		pij_t 	= oldGrid.Cell(it.Top());
		pij_l   = oldGrid.Cell(it.Left());
		pij_r 	= oldGrid.Cell(it.Right());

		A    	= (pij_l+pij_r)/(dx*dx);
		B    	= (pij_d+pij_t)/(dy*dy);

		corr = A+B-rhs->Cell(it);

		grid->Cell(it) = 0.25*(pij_r + pij_l + pij_t + pij_d - (dx*dx)*rhs->Cell(it));


		// Save new p_ij
		// Calculate and summ residual
		real_t lRes = localRes(it,grid,rhs);
		res += lRes;

		// Next cell
		it.Next();
	}

	// Norm residual
	return res/n;
}
