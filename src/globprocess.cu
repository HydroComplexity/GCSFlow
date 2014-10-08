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


// This file provides global functions for general processing for both 2D and 3D
// flow model. Some functions are used to replace built-in CUDA functions. Some
// are used for debugging purposes only.

#include <stdio.h>
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// CopyVariable()
//    Copies data between two variables in the same GPU Device.
//    Input: array var_in
//    Output: array var_out
// --------------------------------------------------------------------
__global__ void CopyVariable(double *var_in, double *var_out, int size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < size) {
    // Transfer data and memory
    var_out[tid] = var_in[tid];

    // Update thread id if vector is long
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// TopForcing()
//    Scale up rainfall to the whole domain
//    Input: ppt (from data)
//    Output: eff_rain (effective rainfall)
// --------------------------------------------------------------------
__global__ void TopForcing(double ppt, double *eff_rain, int size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < size) {
    eff_rain[tid] = ppt;
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// PondHeadInit()
//    Initialize Ponding Head at all points
//    Input: minimum Psi, Pond head (PH)
//    Output: Pond head (PH)
// --------------------------------------------------------------------
__global__ void PondHeadInit(double *ph, int size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < size) {
    ph[tid] = psi_min;
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// TopBottomBound2D()
//    Initialize Ponding Head at all points
//    Input: minimum Psi, Pond head (PH)
//    Output: Pond head (PH)
// --------------------------------------------------------------------
__global__ void TopBottomBound2D(double *Hs, double *Ztopo, double *K2n,
    double *K2s, int BC2D, int M, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    // no-flow BCs
    if (BC2D == 0) {
      Hs[tid] = Hs[N+tid];
      Hs[(M-1)*N+tid] = Hs[(M-2)*N+tid];

    } else {    // Critical depth flow BCs
      Hs[tid] = hcri + Ztopo[tid];
      Hs[(M-1)*N+tid] = hcri + Ztopo[(M-1)*N+tid];
    }

    K2s[tid] = K2s[N+tid];
    K2n[(M-1)*N+tid] = K2n[(M-2)*N+tid];

    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// LeftRightBound2D()
//    Initialize Ponding Head at all points
//    Input: minimum Psi, Pond head (PH)
//    Output: Pond head (PH)
// --------------------------------------------------------------------
__global__ void LeftRightBound2D(double *Hs, double *Ztopo, double *K2e,
    double *K2w, int BC2D, int M, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < M) {
    // no-flow BCs
    if (BC2D == 0) {
      Hs[tid*N] = Hs[tid*N+1];
      Hs[(tid+1)*N-1] = Hs[(tid+1)*N-2];

    } else {    // Critical depth flow BCs
      Hs[tid*N] = hcri + Ztopo[tid*N];
      Hs[(tid+1)*N-1] = hcri + Ztopo[(tid+1)*N-1];
    }

    K2w[tid*N] = K2w[tid*N+1];
    K2e[(tid+1)*N-1] = K2e[(tid+1)*N-2];
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// EstimateVelocity()
//    Estimate 2D velocity on surface in overland flow model
// --------------------------------------------------------------------
__global__ void EstimateVelocity(double *u, double *v, double *K2w, double *K2e,
    double *K2n, double *K2s, double *Hs, double *h, int M, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < M*N) {
    if (h[tid] > hmin) {
      u[tid] = -0.5*(K2w[tid] + K2e[tid])/h[tid] * (Hs[tid+1]-Hs[tid-1])/(2*dx);
      v[tid] = 0.5*(K2n[tid] + K2s[tid])/h[tid] * (Hs[tid+N]-Hs[tid-N])/(2*dy);

    } else {
      u[tid] = 0;
      v[tid] = 0;
    }
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// GetOutlet()
//    Get the depth and discharge at the outlet designed for Benchmarks
//    Input: h, u, v, M, N, P
//    Output: hout, uout, vout
// --------------------------------------------------------------------
__global__ void GetOutlet(double *h, double *houtlet, double *u, double *uout,
    double *v, double *vout, int M, int N, int t) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int ind = 2;
  while (tid < M) {
    houtlet[t*M+tid] = h[(tid+1)*N-ind];
    vout[t*M+tid] = v[(tid+1)*N-ind];
    uout[t*M+tid] = u[(tid+1)*N-ind];
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// getqss()
//    Get the subsurface flow after water balance with Infiltration
//    Input: IN, N, t
//    Output: qss
// --------------------------------------------------------------------
__global__ void getqss(double *IN, double *qss, int N, int t) {
  int I, i, j;
  i = 10; j = 10;
  I = j*N + i;
  qss[t] = IN[I];
}


// --------------------------------------------------------------------
// VarPrint()
//    Print out the values of 2D and 3D variables into sreen for debug
//    For printing 2D variable, use P = 1.
//    Input: Var, M, N, P.
//    Output: 
// --------------------------------------------------------------------
__global__ void VarPrint(double *Var, int M, int N, int P){
  for (int k=0; k < P; k++) {
    for (int i=0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        printf("%4.3f ", Var[k*M*N+i*M+j]);
      }
      printf("\n");
    }
    printf("\n"); printf("\n");
  }
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------