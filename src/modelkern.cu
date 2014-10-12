// Copyright (C) 2014, HydroComplexity Group
// All rights reserved.
//
// GPU-BASED CONJUNCTIVE SURFACE-SUBSURFACE FLOW MODEL (GCSFlow)
// GCSFlow model is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// GCSFlow is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GCSFlow; if not, see <http://www.gnu.org/licenses/>.
//
// Author: levuvietphong@gmail.com (Phong Le)


// This file provides the main CUDA function for GCSFlow model. The functions
// are linked to global and device functions for computation.

#include <string.h>
#include <netcdf.h>     // For using NetCDF format
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/timing.h"
#include "../include/numerical.h"
#include "../include/meminitialize.h"
#include "../include/memtransfer.h"
#include "../include/globprocess.h"
#include "../include/ressave.h"
#include "../include/dataload.h"
#include "../include/class.h"
#include "../include/variableclass.h"
#include "../include/constants.h"
#include "../src/declare.cuh"


// --------------------------------------------------------------------
// NumericalModel()
//    Main code for computation in GPU device. The function includes both C++ 
//    and CUDA commands. NumericalModel function does the following:
//      1. Declare device variables;
//      2. Allocate memory for device variables;
//      3. Transfer data from Host (config file) to device under
//         __constant__ memory;
//      4. Transfer data from Host (load from input file) to device
//         global/shared memory;
//      5. Initialize variables in device memory for initial conditions;
//      6. Run computation per time step - WorkPerStep();
//      7. Save data to NetCDF file and print out information;
// --------------------------------------------------------------------
int NumericalModel(ProjectClass *Project, FileNameClass *File,
    DomainClass *Topo, SoilPropertiesClass *Soil, SubsurfaceClass *Richard3D,
    OverlandClass *Swe2D, HostData2D *Host2D, HostData3D *Host3D,
    HostDataTime *HDTime, HostVarTime *HVTime) {
  int tt, time_data, time_steps, substeps, printstep, savestep;
  int M, N, P, BSZ, TSZ;

  M = Topo->My;     // Dimension in Y direction
  N = Topo->Nx;     // Dimension in X direction
  P = Topo->Pz;     // Dimension in Z direction
  
  // M0, N0. P0 are pre-defined in header files
  if (M!=M0 || N!=N0 || P!=P0){
    printf("The following criteria must be matched: \n");
    printf("M = M0, N = N0, P = P0\n");
    printf("M = %d, M0 = %d, N = %d, N0 = %d, P = %d, P0 = %d\n", 
            M, M0, N, N0, P, P0);
    printf("Change value of M0, N0, P0 in /include/constants.h to match.\n");
    exit(1);
  }
  
  // Assign forcing and options
  time_data = Project->time_data;
  time_steps = Project->time_steps;
  substeps = Project->substeps;
  printstep = Project->PrintSubStep;
  savestep = Project->SaveStep;
  char logfile[32];

  // GPU configuration
  BSZ = Project->BSZ;
  TSZ = Project->TSZ;
  dim3 dimBlockZ(Project->BLOCK_DIMX, Project->BLOCK_DIMY);
  dim3 dimBlockX(Project->BLOCK_DIMY, Project->BLOCK_DIMZ);
  dim3 dimBlockY(Project->BLOCK_DIMX, Project->BLOCK_DIMZ);
  dim3 dimGrid(BSZ, TSZ);
  

  // Allocate CUDA variables in Device
  cudaInitializeData(M, N, P, time_data, time_steps);

  // Transfer data from Host to __constant__ variables in Device
  cudaConstTransfer(Project, Topo, Soil, Richard3D, Swe2D);

  // Transfer data from Host to Device memory
  TransferData(M, N, P, time_data, time_steps, substeps, Host2D, Host3D, HDTime, HVTime);

  // Initialize and copy memory inside GPU
  InitializeVariables(M, N, P, BSZ, TSZ);

  cudaDeviceSynchronize();

  // Measuring time and save to logfile
  double start = GetTime();
  double stop;
  snprintf(logfile, sizeof(char) * 32, "log_%s.txt", Project->Name);
  std::ofstream print_time(logfile, std::ios::app);

  
  // --- TIME LOOP START: --------------------------------------------------
  // The time loop will do the following task sequentially:
  //    1)  Use boundary switching approach or compare infiltrability and 
  //        rainwater supply
  //    2)  Set appropriate boundary condition for Groundwater flow model
  //    3)  Run ADI for 3D Groundwater flow model
  //    4)  Set boundary condition for SWE model
  //    5)  Run ADI/SWE for 2D surface flow model
  // -----------------------------------------------------------------------
  tt = 0; 
  for (int t_data = 0; t_data < time_data; t_data++) {
    for (int it = 0; it < substeps; it++) {
      SimulationPerStep(Project, Swe2D, HVTime, M, N, P, tt, BSZ, TSZ, dimGrid, 
          dimBlockX, dimBlockY, dimBlockZ);      
      
      // Obtain timing per printstep iteration
      if (tt % printstep == 0) {
        stop = GetTime();
        printf("\t Simulation %d of %d completed. Time = %0.3f (s)\n", tt, 
            time_steps, stop - start);
        print_time<< "\t Simulation "<< tt << " of " << time_steps << 
            " completed. Time =" << stop - start << std::endl;
      }
      tt++;
    }

    // Calculate velocity and discharge
    EstimateVelocity<<<TSZ, BSZ>>>(u2_d, v2_d, K2w_d, K2e_d, K2n_d, K2s_d, 
        Hs_out_d, h_d, M, N);

    // Get info at the outlet
    GetOutlet<<<TSZ, BSZ>>>(h_d, hout_d, u2_d, uout_d, v2_d, vout_d, M, N, t_data);

    if (t_data % savestep == 0) {
      //SavePerStep(Project, Host2D, Host3D, File, M, N, P, t_data);
    }
  }

  // --- END OF TIME LOOPS -------------------------------------------------

  print_time.close();
  printf("\t FINISHED...!!!\n\n");
  printf(" PROGRAM END!\n\n\n");

  // Save output at the outlet
  SaveOutlet(File, Swe2D, HVTime, M, N, P, time_data, time_steps);

  // Free GPU memory
  FreeGPUMemory();
  return 0;
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------