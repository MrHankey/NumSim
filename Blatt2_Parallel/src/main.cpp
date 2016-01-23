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

//------------------------------------------------------------------------------
#include "typedef.hpp"
//#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "oclmanager.hpp"
#ifdef USE_DEBUG_VISU
	#include "visu.hpp"
#endif
#include "vtk.hpp"

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>

int main(int argc, char **argv) {

  // Create parameter and geometry instances with default values
  Parameter param;
  Geometry geom;
  OCLManager oclmanager(&geom);

  if ( argc >= 2)
    param.Load(argv[1]);
  else
    param.Load("parameter.txt");

  if ( argc >= 3)
    geom.Load(argv[2]);
  else
    geom.Load("geometry.txt");

  // Create the fluid solver
  Compute comp(&geom, &param, &oclmanager);

  geom.PrintVariables();
  param.PrintVariables();
  // check if folder "VTK" exists
  struct stat info;

  if (stat("VTK", &info) != 0) {
    if( system("mkdir VTK") == -1 )
    {
    	std::cout << "Error creating VTK folder" << std::endl;
    }
  }

// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
  multi_index_t offset_zero_int = {0, 0};
  multi_index_t offset_one_int = {1, 1};
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(800, 800, 1, offset_zero_int, offset_one_int);
#endif // USE_DEBUG_VISU

  // Create a VTK generator;
  // use offset as the domain shift
  multi_real_t offset;
  multi_index_t dims = {1,1};
  offset[0] = 0 * (geom.Mesh()[0] * (real_t)(geom.Size()[0] - 2));
  offset[1] = 0 * (geom.Mesh()[1] * (real_t)(geom.Size()[1] - 2));
  VTK vtk(geom.Mesh(), geom.Size());

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

  oclmanager.Initialize();
  oclmanager.InitFields();

  clock_t begin_full;
  begin_full = clock();
  clock_t begin_step;

  static real_t timing_step = 0.0f;

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
#ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
    case -1:
      return -1;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    case 4:
	  visugrid = comp.GetVorticity();
	  break;
    case 5:
	  visugrid = comp.GetStream();
	  break;
    case 6:
      visugrid = comp.GetT();
      break;
    case 7:
      visugrid = comp.GetRes();
      break;
    default:
      break;
    };
#endif // USE_DEBUG_VISU

    // Create VTK Files in the folder VTK
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
#ifdef VTK_OUTPUT
    vtk.Init("VTK/field");
    //vtk.AddRank();
    vtk.AddField("Cell Velocity", comp.GetU(), comp.GetV());
    //vtk.SwitchToPointData();
    vtk.AddField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddScalar("Pressure", comp.GetP());
    vtk.AddScalar("Vorticity", comp.GetVorticity());
    vtk.AddScalar("Stream", comp.GetStream());
    vtk.AddScalar("Temperature", comp.GetT());
    //vtk.AddPointScalar("Residual", comp.GetRes());
    vtk.Finish();
#endif

    begin_step = clock();
    // Run a few steps
    real_t startTime = comp.GetTime();
    for (uint32_t i = 0; i < 99; ++i)
    //while ( comp.GetTime() - startTime < 0.04)
      comp.TimeStep(false);
    comp.TimeStep(true);
    clock_t end = clock();
	double elapsed_secs = double(end - begin_step) / CLOCKS_PER_SEC;
	timing_step += elapsed_secs;

    //std::cin.ignore();
  }


  clock_t end = clock();
  double elapsed_secs = double(end - begin_full) / CLOCKS_PER_SEC;
  std::cout << "duration: " << elapsed_secs << "step: " << timing_step << " solver_time: " << comp._solver_time << std::endl;
  std::cout << "buf: " << comp._time_buf << " kernel: " << comp._time_kernel << " buf_read: " << comp._time_buf_read << " res: " << comp._time_res << std::endl;


  return 0;
}
