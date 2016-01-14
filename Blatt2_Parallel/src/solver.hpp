/*
 * Copyright (C) 2015   Malte Brunn
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

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __SOLVER_HPP
#define __SOLVER_HPP
//------------------------------------------------------------------------------

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

/** abstract base class for an iterative solver
 */
class Solver {
public:
  /// Constructor of the abstract Solver class
  Solver(const Geometry *geom);
  /// Destructor of the Solver Class
  virtual ~Solver();

  /// This function must be implemented in a child class
  // @param [in][out] grid current values
  // @param [in]      rhs  right hand side values
  // @returns accumulated residual
  virtual real_t Cycle(Grid *grid, const Grid *rhs) const = 0;

protected:
  const Geometry *_geom;

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const;
};

class Jacobi : public Solver {
public:
	/// Constructs an actual Jacobi solver
	  Jacobi(const Geometry *geom);
	  /// Destructor
	  ~Jacobi();

	  /// Returns the total residual and executes a solver cycle
	  // @param grid current pressure values
	  // @param rhs right hand side
	  real_t Cycle(Grid *grid, const Grid *rhs) const;
};

//------------------------------------------------------------------------------
/** concrete SOR solver
 */
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry *geom, const real_t &omega);
  /// Destructor
  ~SOR();

  /// Returns the total residual and executes a solver cycle
  // @param grid current pressure values
  // @param rhs right hand side
  real_t Cycle(Grid *grid, const Grid *rhs) const;

protected:
  real_t _omega;
};
//------------------------------------------------------------------------------

/** concrete Jacobi OpenCL solver
 */
class JacobiOCL{
public:
  /// Constructs an actual SOR solver
  JacobiOCL( const Geometry *geom);
  /// Destructor
  ~JacobiOCL();

  /// Returns the total residual and executes a solver cycle
  // @param grid current pressure values
  // @param rhs right hand side
  real_t Cycle(Grid *grid, const  Grid *rhs);
  void InitializeBuffers();
  void UpdateBuffers(Grid *grid, const Grid *rhs);
  //real_t Cycle(Grid *grid, const Grid *rhs) const;

  double _time_res;
  double _time_kernel;
  double _time_buffer;
  double _time_buffer_read;

protected:
  const Geometry* _geom;
  cl::Kernel _kernel;
  cl::CommandQueue _queue;
  cl::Buffer _bufOld;
  cl::Buffer _bufRHS;
  cl::Buffer _bufNew;
  cl::Buffer _bufLocalResiduals;
  cl::Context _context;
  std::vector<cl::Device> _all_devices;
  cl::Program _program;
};
//------------------------------------------------------------------------------

/** concrete Red or Balck SOR solver
 */
class RedOrBlackSOR : public SOR {
public:
  /// Constructs an actual SOR solver
  RedOrBlackSOR(const Geometry *geom, const real_t &omega);
  /// Destructor
  ~RedOrBlackSOR();

  real_t RedCycle(Grid *grid, const Grid *rhs) const;
  real_t BlackCycle(Grid *grid, const Grid *rhs) const;
};
//------------------------------------------------------------------------------
#endif // __SOLVER_HPP
