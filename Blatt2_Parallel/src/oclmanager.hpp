#include "typedef.hpp"
#include "geometry.hpp"
#include <iostream>
#include <fstream>

#ifndef __OCLMANAGER_HPP
#define __OCLMANAGER_HPP

//#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

class OCLManager {
public:
	OCLManager(Geometry* geom);
	~OCLManager();

	void Initialize();
	void InitFields();
	void SetP(Grid* p);

	real_t ReduceMaxVelocity();
	real_t ReduceResidual();
	//real_t ReduceSum(Grid* grid, index_t length);

	//Kernels
	cl::Kernel _kernel_solver;
	cl::Kernel _kernel_rhs;
	cl::Kernel _kernel_momentumeq;
	cl::Kernel _kernel_newvel;
	cl::Kernel _kernel_reduce_sum;
	cl::Kernel _kernel_reduce_max;

	//Grid Buffers
	cl::Buffer _u;
	cl::Buffer _v;
	cl::Buffer _F;
	cl::Buffer _G;
	cl::Buffer _p;
	cl::Buffer _rhs;
	cl::Buffer _locRes;

	//Constants
	cl::Buffer _h;
	cl::Buffer _h_inv;
	cl::Buffer _hs;
	cl::Buffer _hs_inv;

	cl::CommandQueue _queue;
	cl::Context _context;

private:
	std::vector<cl::Device> _all_devices;
	cl::Program _program;
	cl::Program _program_red;

	Geometry* _geom;

};

inline void checkErr(cl_int err, const char * name) {
    if (err != CL_SUCCESS) {
      std::cerr << "ERROR: " << name  << " (" << err << ")" << std::endl;
      exit(EXIT_FAILURE);
   }
}

#endif
