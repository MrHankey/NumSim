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
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>

int main(int argc, char **argv) {

  // Create parameter and geometry instances with default values
  Communicator comm(&argc, &argv);
  Parameter param;
  Geometry geom(&comm);

  if ( argc >= 2)
    param.Load(argv[1]);
  else
    param.Load("parameter.txt");

  if ( argc >= 3)
    geom.Load(argv[2]);
  else
    geom.Load("geometry.txt");

  // Create the fluid solver
  Compute comp(&geom, &param, &comm);

  if (comm.getRank() == 0) {

    geom.PrintVariables();
	param.PrintVariables();
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }

  int horizontalRes = 800;
  int vertRes = horizontalRes/geom.Length()[0];

// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(horizontalRes / comm.ThreadDim()[0], vertRes / comm.ThreadDim()[1],
              comm.getRank() + 1, comm.ThreadIdx(), comm.ThreadDim());
#endif // USE_DEBUG_VISU

  // Create a VTK generator;
  // use offset as the domain shift
  multi_real_t offset;
  offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0] - 2));
  offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1] - 2));
  VTK vtk(geom.Mesh(), geom.Size(), geom.TotalSize(), offset, comm.getRank(),
          comm.getSize(), comm.ThreadDim());

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

  clock_t begin;
  if ( comm.getRank() == 0)
	  begin = clock();

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
	  visugrid = comp.particelTracingGrid;
	  break;
    case 7:
	  visugrid = comp.StreakLinesGrid;
	  break;
    default:
      break;
    };
#endif // USE_DEBUG_VISU

    // Create VTK Files in the folder VTK
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    vtk.AddPointScalar("Vorticity", comp.GetVorticity());
    vtk.AddPointScalar("Stream", comp.GetStream());
    vtk.AddPointScalar("Tracing", comp.particelTracingGrid);
    vtk.AddPointScalar("StreakLines", comp.StreakLinesGrid);
    vtk.Finish();

    // Run a few steps
    for (uint32_t i = 0; i < 9; ++i)
      comp.TimeStep(false);
    bool printOnlyOnMaster = !comm.getRank();
    comp.TimeStep(printOnlyOnMaster);

    //std::cin.ignore();
  }

  if ( comm.getRank() == 0)
  {
	  clock_t end = clock();
	  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	  std::cout << "duration: " << elapsed_secs << std::endl;
  }

  return 0;
}
