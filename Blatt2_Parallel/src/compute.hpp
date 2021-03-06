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
#include "solver.hpp"

#include <list>
#include <vector>
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(const Geometry *geom, const Parameter *param,
          const Communicator *comm = 0);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  void TimeStep(bool printInfo);

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetRHS() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();

  std::list<multi_real_t> particelTracing;
  std::list<multi_real_t> StreakLines;
  std::vector<multi_real_t> particelTracing_bot;
  std::vector<multi_real_t> StreakLines_bot;

  Grid *particelTracingGrid;
  Grid *StreakLinesGrid;

private:
  // current timestep
  real_t _t;

  // donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  // limit for residual
  real_t _epslimit;

  // velocities
  Grid *_u;
  Grid *_v;

  // pressure
  Grid *_p;

  // temporäre Grids zum debuggen inhalt wird nie benutzt!
  Grid *_temp1;
  Grid *_temp2;

  //Pressure on Right Side and Left Side
  real_t _pL;
  real_t _pR;

  // prel. vel
  Grid *_F;
  Grid *_G;

  // right-hand side
  Grid *_rhs;

  // container for interpolating whichever values
  Grid *_tmp;
  Grid *_tmpVorticity;
  Grid *_tmpStream;

  SOR *_solver;

  const Geometry *_geom;
  const Parameter *_param;
  const Communicator *_comm;

  /// Compute the new velocites u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
