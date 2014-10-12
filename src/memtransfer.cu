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


// This file provides functions to transfer all memory from host to device for
// the simulation.

#include <stdio.h>
#include "../include/class.h"
#include "../include/variableclass.h"
#include "../include/constantextern.h"

// --------------------------------------------------------------------
// OneDprint()
//    Function for printing out one dimensional variable for debugging
// --------------------------------------------------------------------
__global__ void OneDprint() {
  printf("dx, dy, dz, dt: %4.2f, %4.2f, %4.2f, %6.5f \n", dx, dy, dz, dt);
  printf("delta, hmin, hcri, K0: %10.8f, %7.6f, %7.6f, %7.6f \n", 
      delta, hmin, hcri, K0);
  printf("alpha, theta_R, theta_S, nv: %4.2f, %4.2f, %4.2f, %4.2f \n", 
      alpha, theta_R, theta_S, nv);
  printf("poros, Ss, psi_min: %4.2f, %10.8f, %7.6f, %7.6f \n", 
      poros, Ss, psi_min);
  printf("stop_tol, maxiterX, maxiterY, maxiterZ: %10.8f, %d, %d, %d \n", 
      stop_tol, maxiterX, maxiterY, maxiterZ);
}


// --------------------------------------------------------------------
// cudaConstTransfer()
//    Copy memory from configuratio file to CUDA constant memory
// --------------------------------------------------------------------
void cudaConstTransfer(ProjectClass *Project, DomainClass *Topo, 
    SoilPropertiesClass *Soil, SubsurfaceClass *Richard3D, 
    OverlandClass *Swe2D) {
  double dt_h = Project->dt_data/Project->substeps;
  double m = 1.0 - 1.0/Soil->n;
  cudaMemcpyToSymbol(dx, &Topo->GridSizeDx, sizeof(double));
  cudaMemcpyToSymbol(dy, &Topo->GridSizeDy, sizeof(double));
  cudaMemcpyToSymbol(dz, &Topo->GridSizeDz, sizeof(double));
  cudaMemcpyToSymbol(dt, &dt_h, sizeof(double));

  cudaMemcpyToSymbol(delta, &Swe2D->Delta, sizeof(double));
  cudaMemcpyToSymbol(hmin, &Swe2D->MinimumDepth, sizeof(double));
  cudaMemcpyToSymbol(hcri, &Swe2D->CriticalDepth, sizeof(double));
  cudaMemcpyToSymbol(K0, &Swe2D->K0, sizeof(double));

  cudaMemcpyToSymbol(alpha, &Soil->Alpha, sizeof(double));
  cudaMemcpyToSymbol(theta_R, &Soil->Theta_R, sizeof(double));
  cudaMemcpyToSymbol(theta_S, &Soil->Theta_S, sizeof(double));
  cudaMemcpyToSymbol(nv, &Soil->n, sizeof(double));
  cudaMemcpyToSymbol(mv, &m, sizeof(double));
  cudaMemcpyToSymbol(poros, &Soil->porosity, sizeof(double));
  cudaMemcpyToSymbol(Ss, &Soil->Ss, sizeof(double));

  cudaMemcpyToSymbol(psi_min, &Richard3D->Psimin, sizeof(double));
  cudaMemcpyToSymbol(stop_tol, &Richard3D->stoptol, sizeof(double));
  cudaMemcpyToSymbol(Psi_Bottom, &Richard3D->PsiBottom, sizeof(double));
  cudaMemcpyToSymbol(Psi_North, &Richard3D->PsiNorth, sizeof(double));
  cudaMemcpyToSymbol(Psi_South, &Richard3D->PsiSouth, sizeof(double));
  cudaMemcpyToSymbol(Psi_East, &Richard3D->PsiEast, sizeof(double));
  cudaMemcpyToSymbol(Psi_West, &Richard3D->PsiWest, sizeof(double));
  cudaMemcpyToSymbol(BoundBottom, &Richard3D->BoundBottom, sizeof(int));
  cudaMemcpyToSymbol(BoundNorth, &Richard3D->BoundBottom, sizeof(int));
  cudaMemcpyToSymbol(BoundSouth, &Richard3D->BoundSouth, sizeof(int));
  cudaMemcpyToSymbol(BoundEast, &Richard3D->BoundEast, sizeof(int));
  cudaMemcpyToSymbol(BoundWest, &Richard3D->BoundWest, sizeof(int));
  cudaMemcpyToSymbol(maxiterX, &Richard3D->MaxIterX, sizeof(int));
  cudaMemcpyToSymbol(maxiterY, &Richard3D->MaxIterY, sizeof(int));
  cudaMemcpyToSymbol(maxiterZ, &Richard3D->MaxIterZ, sizeof(int));
}


// --------------------------------------------------------------------
// TransferData()
//    Transfer memory from configuratio file to main CUDA memory
// --------------------------------------------------------------------
void TransferData(int M, int N, int P, int time_data, int time_steps,
    int substeps, HostData2D *Host2D, HostData3D *Host3D, HostDataTime *HDTime, 
    HostVarTime *HVTime) {
  double SIZE3D = M*N*P*sizeof(double);
  double SIZE2D = M*N*sizeof(double);
  double SIZE_TIME = time_steps*sizeof(double);

  memcpy(&Host2D->Hs[0], Host2D->Ztopo, SIZE2D);
  int ind = 0;
  for (int i = 0; i < time_data; i++) {
    for (int j = 0; j < substeps; j++) {
      HVTime->PPT[ind] = HDTime->PPT_data[i]/substeps;
      HVTime->ET[ind] = HDTime->ET_data[i]/substeps;
      ind++;
    }
  }

  // --- Copy memory from Host to Device -----------------------
  cudaMemcpy(Psi_in_d, Host3D->Psi_in, SIZE3D, cudaMemcpyHostToDevice);
  cudaMemcpy(Ksat_d, Host3D->Ksat, SIZE3D, cudaMemcpyHostToDevice);
  cudaMemcpy(Hs_in_d, Host2D->Hs, SIZE2D, cudaMemcpyHostToDevice);
  cudaMemcpy(Ztopo_d, Host2D->Ztopo, SIZE2D, cudaMemcpyHostToDevice);
  cudaMemcpy(mann_d, Host2D->mann, SIZE2D, cudaMemcpyHostToDevice);
  cudaMemcpy(PPT_d, HVTime->PPT, SIZE_TIME, cudaMemcpyHostToDevice);
  cudaMemcpy(ET_d, HVTime->ET, SIZE_TIME, cudaMemcpyHostToDevice);
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------