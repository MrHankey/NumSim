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
#include <iostream>

using namespace cl;
using namespace std;

inline void checkErr(cl_int err, const char * name) {
    if (err != CL_SUCCESS) {
      std::cerr << "ERROR: " << name  << " (" << err << ")" << std::endl;
      exit(EXIT_FAILURE);
   }
}

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

/* JacobiOCL solver */
/// Concrete JacobiOCL solver
//  @param geom   get geometry
//  @param omega  get scaling factor
JacobiOCL::JacobiOCL( const Geometry* geom) {

	cl_int err;

	_time_buffer = 0.0;
	_time_buffer_read = 0.0;
	_time_kernel = 0.0;
	_time_kernel = 0.0;

	_geom  = geom;

	vector<Platform> platforms;
	Platform::get(&platforms);

	// Select the default platform and create a context using this platform and the GPU
	cl_context_properties cps[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)(platforms[0])(),
		0
	};
	_context = Context( CL_DEVICE_TYPE_ALL, cps, nullptr, nullptr, &err);
	checkErr(err, "Context::Context()");

	// Get a list of devices on this platform
	//vector<Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	checkErr(platforms[0].getDevices(CL_DEVICE_TYPE_ALL, &_all_devices), "Platform::getDevices()");
	if(_all_devices.size()==0){
		std::cout<<" No devices found. Check OpenCL installation!\n";
		exit(1);
	}

	// Create a command queue and use the first device
	_queue = CommandQueue(_context, _all_devices[0], 0, &err);
	checkErr(err, "CommandQueue::CommandQueue()");

	cout << "Using device: " << _all_devices[0].getInfo<CL_DEVICE_NAME>() << endl;

	// Read source file
	std::ifstream sourceFile("poisson.cl");
	std::string sourceCode(
		std::istreambuf_iterator<char>(sourceFile),
		(std::istreambuf_iterator<char>()));
	Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	// Make program of the source code in the context
	_program = Program(_context, source, &err);
	checkErr(err, "Program::Program()");

	// Build program for these specific devices
	checkErr(_program.build(_all_devices), "Program::build()");
	std::cout<<" Error building: "<<_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_all_devices[0])<<"\n";

	_kernel = Kernel(_program, "poisson_jacobi", &err);
	checkErr(err, "Kernel::Kernel()");

	InitializeBuffers();


}

/// Destructor
JacobiOCL::~JacobiOCL() {}

void JacobiOCL::InitializeBuffers() {
	cl_int err;
	// Create memory buffers
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
#ifdef PINNED_MEMORY
	_bufOld = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, gridSize * sizeof(real_t), nullptr);
	_bufRHS = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, gridSize * sizeof(real_t), nullptr);
	_bufNew = Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, gridSize * sizeof(real_t), nullptr);
	_bufLocalResiduals = Buffer(_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, gridSize * sizeof(real_t), nullptr);
#else
	_bufOld = Buffer(_context, CL_MEM_READ_ONLY, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() old");
	_bufRHS = Buffer(_context, CL_MEM_READ_ONLY, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() rhs");
	_bufNew = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() new");
	_bufLocalResiduals = Buffer(_context, CL_MEM_WRITE_ONLY, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() res");
#endif

}

void JacobiOCL::UpdateBuffers(Grid* grid, const Grid* rhs) {
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

#ifdef PINNED_MEMORY
	// Get exclusive access to the buffer
	real_t* mappedOldGrid = 	(real_t*) _queue.enqueueMapBuffer(_bufOld, CL_TRUE, CL_MAP_READ, 0, gridSize*sizeof(real_t));
	real_t* mappedRHS = 		(real_t*) _queue.enqueueMapBuffer(_bufRHS, CL_TRUE, CL_MAP_READ, 0, gridSize*sizeof(real_t));
	real_t* mappedNewGrid = 	(float*) _queue.enqueueMapBuffer(_bufNew, CL_TRUE, CL_MAP_READ, 0, gridSize*sizeof(real_t));
	//real_t* mappedLocRes = 	(float*) _queue.enqueueMapBuffer(_bufLocalResiduals, CL_TRUE, CL_MAP_READ, 0, gridSize*sizeof(real_t));

    for ( int i = 0; i < gridSize; i++)
    {
    	mappedOldGrid[i] = grid->_data[i];
    	mappedRHS[i] = rhs->_data[i];
    	mappedNewGrid[i] = grid->_data[i];
    }

	// release the buffer back to OpenCL’s control
	_queue.enqueueUnmapMemObject(_bufOld, mappedOldGrid);
	_queue.enqueueUnmapMemObject(_bufRHS, mappedRHS);
	_queue.enqueueUnmapMemObject(_bufNew, mappedNewGrid);
#else
	checkErr(_queue.enqueueWriteBuffer(_bufOld, CL_TRUE, 0, gridSize*sizeof(real_t), grid->_data), "Queue::enqueueWriteBuffer() old");
	checkErr(_queue.enqueueWriteBuffer(_bufRHS, CL_TRUE, 0, gridSize*sizeof(real_t), rhs->_data), "Queue::enqueueWriteBuffer() rhs");
	checkErr(_queue.enqueueWriteBuffer(_bufNew, CL_TRUE, 0, gridSize*sizeof(real_t), grid->_data), "Queue::enqueueWriteBuffer() new");
#endif

}


/// Returns the total residual and executes a solver cycle of SOR
// @param  grid   get grid values
// @param  rhs    get right hand side of equation
// @return Cycle  calculate all new p_ij for one cycle
real_t JacobiOCL::Cycle(Grid* grid, const Grid* rhs) {

	real_t h_square   = _geom->Mesh()[0]*_geom->Mesh()[0];
	real_t h_square_inv = 1.0/h_square;
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	clock_t begin;
	clock_t end;


	Grid localResiduals = Grid(_geom, 0.0f);
	index_t strideSize = 1;

	begin = clock();

	UpdateBuffers(grid, rhs);
	Buffer clHSquare = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square);
	Buffer clHSquareInv = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square_inv);

	_queue.finish();
	end = clock();
	double elapsed_secs_buf = double(end - begin) / CLOCKS_PER_SEC;
	_time_buffer += elapsed_secs_buf;

	begin = clock();


#ifdef __CL_ENABLE_EXCEPTIONS
	try	{
#endif

		// Copy lists A and B to the memory buffers
		//_queue.enqueueWriteBuffer(_bufOld, CL_TRUE, 0, (index_t)gridSize * sizeof(real_t), grid->_data);
		//_queue.enqueueWriteBuffer(_bufRHS, CL_TRUE, 0, (index_t)gridSize * sizeof(real_t), rhs->_data);
		//_queue.enqueueWriteBuffer(_bufNew, CL_TRUE, 0, (index_t)gridSize * sizeof(real_t), grid->_data);

		// Set arguments to kernel
		checkErr(_kernel.setArg(0, _bufOld), "setArg0");
		checkErr(_kernel.setArg(1, _bufRHS), "setArg1");
		checkErr(_kernel.setArg(2, _bufNew), "setArg2");
		//_kernel.setArg(3, sizeof(clDx), &clDx);
		checkErr(_kernel.setArg(3, _bufLocalResiduals), "setArg3");
		checkErr(_kernel.setArg(4, clHSquare), "setArg4");
		checkErr(_kernel.setArg(5, clHSquareInv), "setArg5");
		checkErr(_kernel.setArg(6, sizeof(real_t)*(strideSize+2)*(strideSize+2), NULL), "setArg6");

		// Run the kernel on specific ND range

		NDRange global((_geom->Size()[0] - 2)/strideSize, (_geom->Size()[1] - 2)/strideSize);
		NDRange local(16,16);
		checkErr(_queue.enqueueNDRangeKernel(_kernel, NullRange, global, local), "enqueueNDRangeKernel");
		//_queue.enqueue

		_queue.finish();

		end = clock();
		double elapsed_secs_kernel = double(end - begin) / CLOCKS_PER_SEC;
		_time_kernel += elapsed_secs_kernel;

		begin = clock();


		/*Grid newGrid = Grid(_geom);
		newGrid.Initialize(0.0);*/
		_queue.enqueueReadBuffer(_bufNew, CL_TRUE, 0, gridSize * sizeof(real_t), grid->_data);
		_queue.enqueueReadBuffer(_bufLocalResiduals, CL_TRUE, 0, gridSize * sizeof(real_t), localResiduals._data);
#ifdef __CL_ENABLE_EXCEPTIONS
	} catch(Error error) {
	   std::cout << "Error initializing OpenCL: " << error.what() << "(" << error.err() << ")" << std::endl;
	   exit(1);
	}
#endif

	_queue.finish();

	end = clock();
	double elapsed_secs_buf_read = double(end - begin) / CLOCKS_PER_SEC;
	_time_buffer_read += elapsed_secs_buf_read;


	begin = clock();

	real_t res = 0;
	for ( index_t i = 0; i < gridSize; i++)
	{
		// sum residual

		res += localResiduals._data[i];

	}

	end = clock();
	double elapsed_secs_res = double(end - begin) / CLOCKS_PER_SEC;
	_time_res += elapsed_secs_res;

	// Norm residual
	return res/gridSize;

	/*for(int i = 0; i < LIST_SIZE; i ++)
		 std::cout << A[i] << " + " << B[i] << " = " << C[i] << std::endl; */
}

SOROCL::SOROCL(const Geometry* geom, const real_t& omega)
{
	_geom = geom;
	_omega = omega;

	cl_int err;

	_time_buffer = 0.0;
	_time_buffer_read = 0.0;
	_time_kernel = 0.0;
	_time_kernel = 0.0;

	_geom  = geom;

	vector<Platform> platforms;
	Platform::get(&platforms);

	// Select the default platform and create a context using this platform and the GPU
	cl_context_properties cps[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)(platforms[0])(),
		0
	};
	_context = Context( CL_DEVICE_TYPE_ALL, cps, nullptr, nullptr, &err);
	checkErr(err, "Context::Context()");

	// Get a list of devices on this platform
	//vector<Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	checkErr(platforms[0].getDevices(CL_DEVICE_TYPE_ALL, &_all_devices), "Platform::getDevices()");
	if(_all_devices.size()==0){
		std::cout<<" No devices found. Check OpenCL installation!\n";
		exit(1);
	}

	// Create a command queue and use the first device
	_queue = CommandQueue(_context, _all_devices[0], 0, &err);
	checkErr(err, "CommandQueue::CommandQueue()");

	cout << "Using device: " << _all_devices[0].getInfo<CL_DEVICE_NAME>() << endl;

	// Read source file
	std::ifstream sourceFile("sor_global.cl");
	std::string sourceCode(
		std::istreambuf_iterator<char>(sourceFile),
		(std::istreambuf_iterator<char>()));
	Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	// Make program of the source code in the context
	_program = Program(_context, source, &err);
	checkErr(err, "Program::Program()");

	// Build program for these specific devices
	err = _program.build(_all_devices);
	std::cout<<" Error building: "<<_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_all_devices[0])<<"\n";
	checkErr(err, "Program::build()");

	_kernel = Kernel(_program, "sor", &err);
	checkErr(err, "Kernel::Kernel()");

	InitializeBuffers();
}

SOROCL::~SOROCL()
{

}

void SOROCL::InitializeBuffers()
{
	cl_int err;
	// Create memory buffers
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	_bufGrid = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() grid");
	_bufRHS = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() rhs");
	_bufLocalResiduals = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() res");
}

void SOROCL::UpdateBuffers(Grid* grid, const Grid* rhs, Grid* zeroGrid)
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	checkErr(_queue.enqueueWriteBuffer(_bufGrid, CL_TRUE, 0, gridSize*sizeof(real_t), grid->_data), "Queue::enqueueWriteBuffer() grid");
	checkErr(_queue.enqueueWriteBuffer(_bufRHS, CL_TRUE, 0, gridSize*sizeof(real_t), rhs->_data), "Queue::enqueueWriteBuffer() rhs");
	checkErr(_queue.enqueueWriteBuffer(_bufLocalResiduals, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid->_data), "Queue::enqueueWriteBuffer() res");

}

real_t SOROCL::Cycle(Grid* grid, const Grid* rhs)
{
	real_t h_square   = _geom->Mesh()[0]*_geom->Mesh()[0];
	real_t h_square_inv = 1.0/h_square;

	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	clock_t begin;
	clock_t end;

	//grid for saving local residuals
	//initialize as 0
	Grid localResiduals = Grid(_geom);
	localResiduals.Initialize(0.0f);

	index_t localSize = 16;
	index_t localCellCount = (localSize+2)*(localSize+2);

	begin = clock();

	//copy grid data to device
	UpdateBuffers(grid, rhs, &localResiduals);

	Buffer clHSquare = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square);
	Buffer clHSquareInv = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square_inv);

	_queue.finish();
	end = clock();
	double elapsed_secs_buf = double(end - begin) / CLOCKS_PER_SEC;
	_time_buffer += elapsed_secs_buf;

	begin = clock();


#ifdef __CL_ENABLE_EXCEPTIONS
	try	{
#endif

		// Set arguments to kernel
		checkErr(_kernel.setArg(0, _bufGrid), "setArg0");
		checkErr(_kernel.setArg(1, _bufRHS), "setArg1");
		checkErr(_kernel.setArg(2, _bufLocalResiduals), "setArg2");
		checkErr(_kernel.setArg(3, clHSquare), "setArg3");
		checkErr(_kernel.setArg(4, clHSquareInv), "setArg4");
		//checkErr(_kernel.setArg(5, sizeof(real_t)*localCellCount, NULL), "setArg5");
		//checkErr(_kernel.setArg(6, sizeof(real_t)*localCellCount, NULL), "setArg6");

		// Run the kernel on specific ND range
		//cout << _geom->Size()[0] - 2 << " " << (_geom->Size()[1] - 2) << endl;
		NDRange global((_geom->Size()[0] - 2), (_geom->Size()[1] - 2));
		NDRange local(localSize,localSize);
		checkErr(_queue.enqueueNDRangeKernel(_kernel, NullRange, global, local), "enqueueNDRangeKernel");

		_queue.finish();

		end = clock();
		double elapsed_secs_kernel = double(end - begin) / CLOCKS_PER_SEC;
		_time_kernel += elapsed_secs_kernel;

		begin = clock();

		_queue.enqueueReadBuffer(_bufGrid, CL_TRUE, 0, gridSize * sizeof(real_t), grid->_data);
		_queue.enqueueReadBuffer(_bufLocalResiduals, CL_TRUE, 0, gridSize * sizeof(real_t), localResiduals._data);
#ifdef __CL_ENABLE_EXCEPTIONS
	} catch(Error error) {
	   std::cout << "Error initializing OpenCL: " << error.what() << "(" << error.err() << ")" << std::endl;
	   exit(1);
	}
#endif

	_queue.finish();

	end = clock();
	double elapsed_secs_buf_read = double(end - begin) / CLOCKS_PER_SEC;
	_time_buffer_read += elapsed_secs_buf_read;


	begin = clock();


	real_t res = 0;
	for ( index_t i = 0; i < gridSize; i++)
	{
		// sum residual

		float val = localResiduals._data[i];
		if ( std::isfinite(val))
		{
			res += val;
		}

	}

	/*InteriorIterator it = InteriorIterator(_geom);
	while (it.Valid())
	{
		real_t loc = localResiduals.Cell(it);
		//if ( loc > 0.0)
			//cout << loc << endl;

		res += loc;
		it.Next();
	}*/

	//cin.ignore();

	end = clock();
	double elapsed_secs_res = double(end - begin) / CLOCKS_PER_SEC;
	_time_res += elapsed_secs_res;

	// Norm residual
	return res/gridSize;
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
