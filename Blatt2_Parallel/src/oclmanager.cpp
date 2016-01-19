#include "oclmanager.hpp"
#include "grid.hpp"

using namespace std;
using namespace cl;

OCLManager::OCLManager(Geometry* geom) {

	_geom = geom;

	cl_int err;

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

	_kernel_solver = Kernel(_program, "sor", &err);
	checkErr(err, "Kernel::Kernel() sor");

	_kernel_newvel = Kernel(_program, "newvel", &err);
	checkErr(err, "Kernel::Kernel() newvel");

	_kernel_rhs = Kernel(_program, "rhs", &err);
	checkErr(err, "Kernel::Kernel() rhs");

	_kernel_momentumeq = Kernel(_program, "momentumeq", &err);
	checkErr(err, "Kernel::Kernel() momentumeq");

	//Initialize();
	//InitFields();
}

OCLManager::~OCLManager() {
}

void OCLManager::InitFields()
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];

	Grid zeroGrid = Grid(_geom);
	zeroGrid.Initialize(0.0f);

	checkErr(_queue.enqueueWriteBuffer(_p, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer p");
	checkErr(_queue.enqueueWriteBuffer(_u, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer u");
	checkErr(_queue.enqueueWriteBuffer(_v, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer v");
	checkErr(_queue.enqueueWriteBuffer(_F, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer F");
	checkErr(_queue.enqueueWriteBuffer(_G, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer G");
	checkErr(_queue.enqueueWriteBuffer(_rhs, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer rhs");
	checkErr(_queue.enqueueWriteBuffer(_locRes, CL_TRUE, 0, gridSize*sizeof(real_t), zeroGrid._data), "manager write buffer locRes");

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

void OCLManager::SetP(Grid* p)
{
	index_t gridSize = _geom->Size()[0]*_geom->Size()[1];
	checkErr(_queue.enqueueWriteBuffer(_p, CL_TRUE, 0, gridSize*sizeof(real_t), p->_data), "ocl amanager set pu");
}
