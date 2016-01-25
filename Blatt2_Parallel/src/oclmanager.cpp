#include "oclmanager.hpp"
#include "grid.hpp"
#include <cmath>

using namespace std;
using namespace cl;

//intel: 15
//k20x: 96
#define NUM_WORKGROUPS 96

OCLManager::OCLManager(Geometry* geom, index_t device) {

	_geom = geom;

	cl_int err;

	std::vector<Platform> platforms;
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
	_queue = CommandQueue(_context, _all_devices[device], 0, &err);
	checkErr(err, "CommandQueue::CommandQueue()");

	cout << "Using device: " << _all_devices[device].getInfo<CL_DEVICE_NAME>() << endl;

	_numWorkGroups = NUM_WORKGROUPS;

	// Read source file
	std::ifstream sourceFile("computations.cl");
	std::string sourceCode(
		std::istreambuf_iterator<char>(sourceFile),
		(std::istreambuf_iterator<char>()));
	Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	// Make program of the source code in the context
	_program = Program(_context, source, &err);
	checkErr(err, "Program::Program()");

	// Build program for these specific devices
	err = _program.build(_all_devices);
	if ( err != 0 )
		std::cout<<" Error building: "<<_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_all_devices[device])<<"\n";

	checkErr(err, "Program::build()");

	_kernel_solver = Kernel(_program, "sor", &err);
	checkErr(err, "Kernel::Kernel() sor");

	_kernel_newvel = Kernel(_program, "newvel", &err);
	checkErr(err, "Kernel::Kernel() newvel");

	_kernel_rhs = Kernel(_program, "rhs", &err);
	checkErr(err, "Kernel::Kernel() rhs");

	_kernel_momentumeq = Kernel(_program, "momentumeq", &err);
	checkErr(err, "Kernel::Kernel() momentumeq");

	_kernel_T = Kernel(_program, "compT", &err);
	checkErr(err, "Kernel::Kernel() T");

	std::ifstream sourceFileRed("reductions.cl");
	std::string sourceCodeRed(
			std::istreambuf_iterator<char>(sourceFileRed),
			(std::istreambuf_iterator<char>()));
	Program::Sources sourceRed(1, std::make_pair(sourceCodeRed.c_str(), sourceCodeRed.length()+1));

	// Make program of the source code in the context
	_program_red = Program(_context, sourceRed, &err);
	checkErr(err, "Program::Program(red)");

	// Build program for these specific devices
	err = _program_red.build(_all_devices);
	if (err != 0)
		std::cout<<" Error building: "<<_program_red.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_all_devices[device])<<"\n";

	checkErr(err, "Program::build(red)");

	_kernel_reduce_sum = Kernel(_program_red, "reduce_sum", &err);
	checkErr(err, "Kernel::Kernel() red_sum");

	_kernel_reduce_max = Kernel(_program_red, "reduce_max", &err);
	checkErr(err, "Kernel::Kernel() red_max");

	//Initialize();
	//InitFields();
}

OCLManager::~OCLManager() {
}

void OCLManager::InitFields()
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	Grid zeroGrid = Grid(_geom);
	zeroGrid.Initialize(0.0001f);

	checkErr(_queue.enqueueWriteBuffer(_p, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer p");
	checkErr(_queue.enqueueWriteBuffer(_u, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer u");
	checkErr(_queue.enqueueWriteBuffer(_v, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer v");
	checkErr(_queue.enqueueWriteBuffer(_F, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer F");
	checkErr(_queue.enqueueWriteBuffer(_G, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer G");
	checkErr(_queue.enqueueWriteBuffer(_rhs, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer rhs");
	checkErr(_queue.enqueueWriteBuffer(_locRes, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer locRes");
	checkErr(_queue.enqueueWriteBuffer(_T, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer T");

}

void OCLManager::Initialize() {
	cl_int err;
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	_u = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() u");
	_v = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() v");
	_p = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() p");
	_F = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() F");
	_G = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() G");
	_rhs = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() G");
	_locRes = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() res");
	_T = Buffer(_context, CL_MEM_READ_WRITE, gridSize * sizeof(real_t), nullptr, &err);
	checkErr(err, "Buffer::Buffer() T");

	real_t h = _geom->Mesh()[0];
	real_t h_inv = 1.0/h;

	real_t h_square   = _geom->Mesh()[0]*_geom->Mesh()[0];
	real_t h_square_inv = 1.0/h_square;

	_hs = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square, &err);
	checkErr(err, "Buffer::Buffer() hs");
	_hs_inv = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_square_inv, &err);
	checkErr(err, "Buffer::Buffer() hsi");

	_h = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h, &err);
	checkErr(err, "Buffer::Buffer() h");
	_h_inv = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(real_t), &h_inv, &err);
	checkErr(err, "Buffer::Buffer() hi");


}

real_t OCLManager::ReduceMaxVelocity()
{
	//intel 15
	cl_uint gridSize1D = _geom->Size()[0] * (_geom->Size()[1] - 1);
	real_t result [NUM_WORKGROUPS];

	NDRange global_red(_numWorkGroups*128);
	NDRange local_red(128);

	Buffer clLength = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_uint), &gridSize1D);
	Buffer clResult = Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(real_t)*_numWorkGroups, &result);

	checkErr(_kernel_reduce_max.setArg(0, _u), "setArg0 red max u");
	checkErr(_kernel_reduce_max.setArg(1, sizeof(real_t)*128, nullptr), "setArg1 red max u");
	checkErr(_kernel_reduce_max.setArg(2, clLength), "setArg2 red max u");
	checkErr(_kernel_reduce_max.setArg(3, clResult), "setArg3 red max u");

	checkErr(_queue.enqueueNDRangeKernel(_kernel_reduce_max, NullRange, global_red, local_red), "enqueueNDRangeKernelSolver");
	_queue.enqueueReadBuffer(clResult, CL_TRUE, 0, sizeof(real_t)*_numWorkGroups, &result);

	real_t uMax = result[0];

	for (index_t i = 1; i < _numWorkGroups; i++)
	{
		uMax = fmax(result[i], uMax);
	}

	checkErr(_kernel_reduce_max.setArg(0, _v), "setArg0 red max v");
	checkErr(_kernel_reduce_max.setArg(1, sizeof(real_t)*128, nullptr), "setArg1 red max v");
	checkErr(_kernel_reduce_max.setArg(2, clLength), "setArg2 red max v");
	checkErr(_kernel_reduce_max.setArg(3, clResult), "setArg3 red max v");

	checkErr(_queue.enqueueNDRangeKernel(_kernel_reduce_max, NullRange, global_red, local_red), "enqueueNDRangeKernelSolver");
	_queue.enqueueReadBuffer(clResult, CL_TRUE, 0, sizeof(real_t)*_numWorkGroups, &result);

	real_t vMax = result[0];
	for (index_t i = 1; i < _numWorkGroups; i++)
	{
		vMax = fmax(result[i], vMax);
	}

	return fmax(uMax,vMax);
}

real_t OCLManager::ReduceResidual()
{

	cl_uint gridSize1D = _geom->Size()[0] * (_geom->Size()[1]);
	real_t result[NUM_WORKGROUPS];


	NDRange global_red(_numWorkGroups*128);
	NDRange local_red(128);

	Buffer clLength = Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_uint), &gridSize1D);
	Buffer clResult = Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(real_t)*_numWorkGroups, &result);

	checkErr(_kernel_reduce_sum.setArg(0, _locRes), "setArg0 red sum");
	checkErr(_kernel_reduce_sum.setArg(1, sizeof(real_t)*128, nullptr), "setArg1 red sum");
	checkErr(_kernel_reduce_sum.setArg(2, clLength), "setArg2 red sum");
	checkErr(_kernel_reduce_sum.setArg(3, clResult), "setArg3 red sum");

	checkErr(_queue.enqueueNDRangeKernel(_kernel_reduce_sum, NullRange, global_red, local_red), "enqueueNDRangeKernelSolver");
	_queue.enqueueReadBuffer(clResult, CL_TRUE, 0, sizeof(real_t)*_numWorkGroups, &result);

	real_t sum = 0.0f;

	for (index_t i = 0; i < _numWorkGroups; i++)
	{
		sum += result[i];
	}

	return sum/gridSize1D;
}

void OCLManager::SetP(Grid* p)
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	checkErr(_queue.enqueueWriteBuffer(_p, CL_TRUE, 0, gridSize*sizeof(real_t), p->_data), "ocl amanager set pu");
}