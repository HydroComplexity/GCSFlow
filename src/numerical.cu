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


// In this file, we provide two functions to save result over the whole domain 
// and at the specified outlet at particular time step chosen by user. The third
// function is developed for computational and numerical solution.

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
// SavePerStep()
//    Save model output over the whole domain per a number of time steps
// --------------------------------------------------------------------
void SavePerStep(ProjectClass *Project, HostData2D *Host2D, HostData3D *Host3D, 
    FileNameClass *File, int M, int N, int P, int t_data) {
  char file2D[32], file3D[32], file_veloc[32], file_IN[32], fdiff[32];
  char iterZ[32], iterX[32], iterY[32];

  if (Project->Save2DOverland == 1) {
    cudaMemcpy(Host2D->h_out, htop, M*N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Host2D->Hs_out, Hs_in_d, M*N*sizeof(double), 
        cudaMemcpyDeviceToHost);
    snprintf(file2D, 32, "%s_%d.nc", File->output2D, t_data);
    SaveOutput2D(file2D, N, M, "h_out", Host2D->h_out, "Hs_out",Host2D->Hs_out);
  }

  if (Project->Save3DSubsurface == 1) {
    cudaMemcpy(Host3D->Psi_out, Psi_out_d, M*N*P*sizeof(double), 
        cudaMemcpyDeviceToHost);
    cudaMemcpy(Host3D->theta_out, theta_out_d, M*N*P*sizeof(double), 
        cudaMemcpyDeviceToHost);
    snprintf(file3D, 32, "%s_%d.nc", File->output3D, t_data);
    SaveOutput3D(file3D, N, M, P, "theta_out", Host3D->theta_out, "Psi_out", 
        Host3D->Psi_out);
  }

  if (Project->SaveIteration == 1) {
    cudaMemcpy(Host2D->iter_z_h, iter_z_d, M*N*sizeof(int), 
        cudaMemcpyDeviceToHost);
    cudaMemcpy(Host2D->iter_x_h, iter_x_d, M*P*sizeof(int), 
        cudaMemcpyDeviceToHost);
    cudaMemcpy(Host2D->iter_y_h, iter_y_d, N*P*sizeof(int), 
        cudaMemcpyDeviceToHost);
    snprintf(iterZ, 32, "%s_%d.nc", "iterZ", t_data);
    snprintf(iterX, 32, "%s_%d.nc", "iterX", t_data);
    snprintf(iterY, 32, "%s_%d.nc", "iterY", t_data);
    SaveOneIntOutput2D(iterZ, N, M, "iterZ", Host2D->iter_z_h);
    SaveOneIntOutput2D(iterX, M, P, "iterX", Host2D->iter_x_h);
    SaveOneIntOutput2D(iterY, N, P, "iterY", Host2D->iter_y_h);
  }

  if (Project->SaveInfiltration == 1){
    cudaMemcpy(Host2D->IN_h, IN, M*N*sizeof(double), cudaMemcpyDeviceToHost);
    snprintf(file_IN, 32, "%s_%d.nc", "Infil", t_data);
    SaveOneOutput2D(file_IN, N, M, "IN", Host2D->IN_h);
  }

  /*
  //cudaMemcpy(Host2D->u2_h, u2_d, M*N*sizeof(double), cudaMemcpyDeviceToHost);
  //cudaMemcpy(Host2D->v2_h, v2_d, M*N*sizeof(double), cudaMemcpyDeviceToHost);
  //cudaMemcpy(Host2D->Psidiff_h, Psidiff, M*N*sizeof(double), cudaMemcpyDeviceToHost);

  snprintf(file_veloc, 32, "%s_%d.nc", File->output_veloc2D, t_data);
  snprintf(fdiff, 32, "%s_%d.nc", "Psidiff", t_data);  

  SaveOutput2D(file_veloc, N, M, "U", Host2D->u2_h, "IN", Host2D->IN_h);  
  SaveOneOutput2D(fdiff, N, M, "Psidiff", Host2D->Psidiff_h);
  */  
}


// --------------------------------------------------------------------
// SaveOutlet()
//    Save model output at the outlet per a number of time steps
// --------------------------------------------------------------------
void SaveOutlet(FileNameClass *File, OverlandClass *Swe2D, HostVarTime *HVTime, 
    int M, int N, int P, int time_data, int time_steps){
  double *hout = new double[M*time_data];
  double *vout = new double[M*time_data];
  double *uout = new double[M*time_data];
  char fhout[32], fuout[32], fvout[32], qss_out[32], numstr[32];

  cudaMemcpy(hout, hout_d, M*time_data*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(vout, vout_d, M*time_data*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(uout, uout_d, M*time_data*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(HVTime->qss, qss_d, time_steps*sizeof(double), 
      cudaMemcpyDeviceToHost);
  
  if (Swe2D->Numeric == 0) {
    snprintf(numstr, 32,"exp");
  } else {
    snprintf(numstr, 32,"adi");
  }
  snprintf(fhout, 32, "%s_%s.nc", File->output_h, numstr);
  snprintf(fuout, 32, "%s_%s.nc", File->output_u, numstr);
  snprintf(fvout, 32, "%s_%s.nc", File->output_v, numstr);
  snprintf(qss_out, 32, "%s.nc", File->output_qss);  

  SaveVariables1D(fhout, M*time_data, "hout", hout);
  SaveVariables1D(fuout, M*time_data, "uout", uout);
  SaveVariables1D(fvout, M*time_data, "vout", vout);
  SaveVariables1D(qss_out, time_steps, "qss", HVTime->qss);
  delete[] hout; delete[] vout; delete[] uout;
}


// --------------------------------------------------------------------
// SimulationPerStep()
//    This function call numerical functions for simulation at every
//    timestep in the model. Simulation includes both 2D overland flow
//    and 3D subsurface flow.
// --------------------------------------------------------------------
void SimulationPerStep(ProjectClass *Project, OverlandClass *Swe2D, 
    HostVarTime *HVTime, int M, int N, int P, int tt, int BSZ, int TSZ,
    dim3 dimGrid, dim3 dimBlockX, dim3 dimBlockY, dim3 dimBlockZ) {
  int BC2D = Project->BC2D;
  double SIZE = M*N*P*sizeof(double);
  // Update flux and rainfall
  TopForcing<<<TSZ, BSZ>>>(HVTime->PPT[tt], eRF, M*N);

  //. . .2D SURFACE FLOW MODELING . . . . . . . . . . . . . . . .
  // Shallow Water Equation - Explicit scheme
  if (Swe2D->Numeric == 0) {
    SWE_Explicit<<<dimGrid, dimBlockZ>>> (Hs_in_d, h_d, Hs_out_d, K2w_d, K2e_d, 
        K2n_d, K2s_d, Ztopo_d, mann_d, eRF, IN, ET_d, M, N, tt);

    cudaMemcpy(Hs_in_d, Hs_out_d, M*N*sizeof(double), cudaMemcpyDeviceToDevice);
    LeftRightBound2D<<<TSZ, BSZ>>>(Hs_in_d, Ztopo_d, K2e_d, K2w_d, BC2D, M, N);
    TopBottomBound2D<<<TSZ, BSZ>>>(Hs_in_d, Ztopo_d, K2n_d, K2s_d, BC2D, M, N);
    cudaDeviceSynchronize();

  } else {  // Shallow Water Equation - Alternating Direction Implicit
    if (tt % 2 == 0) {
      // If time is even, do X-direction first, then Y-direction
      SweX<<<TSZ, BSZ>>>(Hs_in_d, h_d, Hs_out_d, K2w_d, K2e_d, Ztopo_d, mann_d, 
          eRF, IN, ET_d, M, N, tt);
      SweY<<<TSZ, BSZ>>>(Hs_out_d, h_d, Hs_in_d, K2n_d, K2s_d, Ztopo_d, mann_d,
          eRF, IN, ET_d, M, N, tt);
    } else {
      // If time is odd, do Y-direction first, then X-direction
      SweY<<<TSZ, BSZ>>>(Hs_in_d, h_d, Hs_out_d, K2n_d, K2s_d, Ztopo_d, mann_d, 
          eRF, IN, ET_d, M, N, tt);
      SweX<<<TSZ, BSZ>>>(Hs_out_d, h_d, Hs_in_d, K2w_d, K2e_d, Ztopo_d, mann_d, 
          eRF, IN, ET_d, M, N, tt);
    }
    cudaDeviceSynchronize();
  }

  //. . .3D SUBSURFACE FLOW MODELING . . . . . . . . . . . . . .
  // Richards Equation - Alternating Direction Implicit
  RichZ<<<dimGrid, dimBlockZ>>> (Psi_in_d, theta_in_d, K_in_d, Psi_out_d, 
      theta_out_d, K_out_d, Ksat_d, h_d, HVTime->PPT[tt], IN, Psidiff, iter_z_d,
      M, N, P);

  if (tt % 2 == 0) {
    // If time is even, do X-direction first, then Y-direction
    RichX<<<dimGrid, dimBlockX>>>(Psi_out_d, theta_out_d, K_out_d, Psi_in_d, 
        theta_in_d, K_in_d, Ksat_d, iter_x_d, M, N, P);     
    RichY<<<dimGrid, dimBlockY>>>(Psi_in_d, theta_in_d, K_in_d, Psi_out_d, 
        theta_out_d, K_out_d, Ksat_d, iter_y_d, M, N, P);
  } else {
    // If time is odd, do Y-direction first, then X-direction
    RichY<<<dimGrid, dimBlockX>>>(Psi_out_d, theta_out_d, K_out_d, Psi_in_d, 
        theta_in_d, K_in_d, Ksat_d, iter_y_d, M, N, P);
    RichX<<<dimGrid, dimBlockY>>>(Psi_in_d, theta_in_d, K_in_d, Psi_out_d, 
        theta_out_d, K_out_d, Ksat_d, iter_x_d, M, N, P);
  }

  // Synchronize all threads here for memory swap
  cudaDeviceSynchronize();    
  cudaMemcpy(Psi_in_d, Psi_out_d, SIZE, cudaMemcpyDeviceToDevice);
  cudaMemcpy(K_in_d, K_out_d, SIZE, cudaMemcpyDeviceToDevice);
  cudaMemcpy(theta_in_d, theta_out_d, SIZE, cudaMemcpyDeviceToDevice);
  getqss<<<1, 1>>>(IN, qss_d, N, tt);
}
