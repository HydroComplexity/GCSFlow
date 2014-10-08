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


// This file provides functions to initialize and free all memory in device for 
// the  simulations.

#include <string.h>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "../include/globoverland.h"
#include "../include/globsubsurface.h"
#include "../include/globprocess.h"
#include "../include/ressave.h"
#include "../include/class.h"
#include "../include/variableclass.h"
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// cudaInitializeData()
//    Initialize main variables in main CUDA memory for both overland flow and
//    subsurface flow models.
// --------------------------------------------------------------------
void cudaInitializeData(int M, int N, int P, int time_data, int time_steps) {
  double SIZE3D = M*N*P*sizeof(double);
  double SIZE2D = M*N*sizeof(double);
  // Subsurface model variables
  cudaMalloc((void**)&Psi_in_d, SIZE3D);
  cudaMalloc((void**)&theta_in_d, SIZE3D);
  cudaMalloc((void**)&K_in_d, SIZE3D);
  cudaMalloc((void**)&C_in_d, SIZE3D);
  cudaMalloc((void**)&Psi_out_d, SIZE3D);
  cudaMalloc((void**)&theta_out_d, SIZE3D);
  cudaMalloc((void**)&K_out_d, SIZE3D);
  cudaMalloc((void**)&Ksat_d, SIZE3D);

  cudaMalloc((void**)&eRF, SIZE2D);
  cudaMalloc((void**)&IN, SIZE2D);
  cudaMalloc((void**)&Psidiff, SIZE2D);
  cudaMalloc((void**)&hbottom, SIZE2D);
  cudaMalloc((void**)&iter_z_d, M*N*sizeof(int));
  cudaMalloc((void**)&iter_x_d, M*P*sizeof(int));
  cudaMalloc((void**)&iter_y_d, N*P*sizeof(int));

  // Overland - Surface model variables
  cudaMalloc((void**)&Hs_in_d, SIZE2D);
  cudaMalloc((void**)&h_d, SIZE2D);
  cudaMalloc((void**)&Hs_out_d, SIZE2D);
  cudaMalloc((void**)&Ztopo_d, SIZE2D);
  cudaMalloc((void**)&mann_d, SIZE2D);
  cudaMalloc((void**)&K2w_d, SIZE2D);
  cudaMalloc((void**)&K2e_d, SIZE2D);
  cudaMalloc((void**)&K2n_d, SIZE2D);
  cudaMalloc((void**)&K2s_d, SIZE2D);
  cudaMalloc((void**)&u2_d, SIZE2D);
  cudaMalloc((void**)&v2_d, SIZE2D);
  
  cudaMalloc((void**)&PPT_d, time_steps*sizeof(double));
  cudaMalloc((void**)&ET_d, time_steps*sizeof(double));
  cudaMalloc((void**)&hout_d, M*time_data*sizeof(double));
  cudaMalloc((void**)&vout_d, M*time_data*sizeof(double));
  cudaMalloc((void**)&uout_d, M*time_data*sizeof(double));
  cudaMalloc((void**)&qss_d, time_steps*sizeof(double));

  // Set variables to 0 for initialization
  cudaMemset(eRF, 0, SIZE2D);
  cudaMemset(IN, 0, SIZE2D);
  cudaMemset(K2w_d, 0, SIZE2D);
  cudaMemset(K2e_d, 0, SIZE2D);
  cudaMemset(K2n_d, 0, SIZE2D);
  cudaMemset(K2s_d, 0, SIZE2D);
  cudaMemset(u2_d, 0, SIZE2D);
  cudaMemset(v2_d, 0, SIZE2D);

  cudaMemset(iter_z_d, 0, M*N*sizeof(int));
  cudaMemset(iter_x_d, 0, M*P*sizeof(int));
  cudaMemset(iter_y_d, 0, N*P*sizeof(int));
}


// --------------------------------------------------------------------
// InitializeVariables()
//    Initialize variables and initial conditions in overland flow and 
//    subsurface flow models
// --------------------------------------------------------------------
void InitializeVariables(int M, int N, int P, int BSZ, int TSZ) {
  double SIZE3D = M*N*P*sizeof(double);
  double SIZE2D = M*N*sizeof(double);
  vanGenuchtenIntial<<<TSZ, BSZ>>>(theta_in_d, K_in_d, Ksat_d, Psi_in_d, M*N*P);
  PondHeadInit<<<TSZ, BSZ>>>(h_d, M*N);
  SweHInit<<<TSZ, BSZ>>>(Hs_in_d, Ztopo_d, h_d, M*N);
  cudaMemcpy(Psi_out_d, Psi_in_d, SIZE3D, cudaMemcpyDeviceToDevice);
  cudaMemcpy(K_out_d, K_in_d, SIZE3D, cudaMemcpyDeviceToDevice);
  cudaMemcpy(theta_out_d, theta_in_d, SIZE3D, cudaMemcpyDeviceToDevice);
  cudaMemcpy(Hs_out_d, Hs_in_d, SIZE2D, cudaMemcpyDeviceToDevice);
}


// --------------------------------------------------------------------
// FreeGPUMemory()
//    Free device memory at the end of the simulation
// --------------------------------------------------------------------
void FreeGPUMemory() {
  cudaFree(Hs_in_d); cudaFree(h_d); cudaFree(Hs_out_d); cudaFree(Ztopo_d);
  cudaFree(Psi_in_d); cudaFree(Psi_out_d); cudaFree(K_in_d); cudaFree(K_out_d);
  cudaFree(C_in_d); cudaFree(vout_d); cudaFree(uout_d); cudaFree(hout_d);
  cudaFree(eRF); cudaFree(PPT_d); cudaFree(ET_d); cudaFree(IN);
  cudaFree(K2w_d); cudaFree(K2e_d); cudaFree(K2n_d); cudaFree(K2s_d);
  cudaFree(u2_d); cudaFree(v2_d); cudaFree(qss_d); cudaFree(hbottom);
}


// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------