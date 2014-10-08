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


// This file provides global functions for overland flow model. Both explicit 
// and Alternating Direction Implicit methods are implemented. Global functions
// are linked with device functions 

#include "../include/constants.h"
#include "../include/devmaxfind.h"
#include "../include/devtdma.h"
#include "../include/devoverland.h"
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// SweHInit()
//    Initialize water depth values in the domain in GPU/device memory
//    Input...
//        * var_in1: input variable 1 (Water elevation)
//        * var_in2: input variable 2 (Ground elevation)
//        * size: size (M*N in 1D) of the domain
//    Output...
//        * var_out: output variable (water depth)
// --------------------------------------------------------------------
__global__ void SweHInit(double *var_in1, double *var_in2, double *var_out,
                         int size) {
  // Get thread id
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < size) {
    // Transfer data and memory and calculation
    var_out[tid] = var_in1[tid] - var_in2[tid];

    // Thread id update
    tid += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// SweX()
// SweY()
//    These functions sets up the discritized equation for solving tridiagonal
//    systems in n 2D overland flow model using ADI method.
//    SweX() for horizontal x direction modeling in ADI.
//    SweY() for horizontal y direction modeling in ADI.
//    References:
//        * Wasantha Lal 1998, J. Hydraulic Engineering
//        * Morita and Yeh (2002) Int. J. Numer. Meth. Fluids
//
//    Input...
//        * Hs_in: Water elevation [L]
//        * Ztopo: Ground surface elevation [L]
//        * h: Water depth [L]
//        * mann: Manning's coefficient [L^{-1/3} T]
//        * eRF: effective rainfall [L/T]
//        * ET: Evaporation [L/T]
//        * Scalar indexes
//
//    Output...
//        Matrix A and vector RHS d
//        * Hs_out: Water elevation[L]
//        * h: Water depth [L]
//        * IN: Infiltration [L/T]
//        * K2n and K2s: K values in north and south directions
//        * K2w and K2e: K values in east and west directions
// --------------------------------------------------------------------
__global__ void SweX(double *Hs_in, double *h, double *Hs_out, double *K2w,
                     double *K2e, double *Ztopo, double *mann, double *eRF,
                     double *IN, double *ET, int M, int N, int t) {
  int j0 = threadIdx.x + blockIdx.x * blockDim.x;
  int j = j0+1;

  // Declaration of local variable in each thread
  double a_x[N0-2], b_x[N0-2], c_x[N0-2], d_x[N0-2], f_x[N0-2];
  double c_star_x[N0-2], d_star_x[N0-2];

  while (j > 0 && j < M-1) {
    // Set the tri-diagonal systems in X-direction
    SweXMatrix(a_x, b_x, c_x, d_x, N, j, t, Hs_in, h, mann, eRF, IN, ET, K2w, K2e);

    // Solve the tri-diagonal systems
    ThomasAlgorithm(a_x, b_x, c_x, d_x, f_x, c_star_x, d_star_x, N-2);

    // Update global values from TDMA solution
    for (int i = 1; i < N-1; i++) {
      Hs_out[j*N+i] = f_x[i-1];
      h[j*N+i] = Hs_out[j*N+i] - Ztopo[j*N+i];
    }

    j += gridDim.x * blockDim.x;
  }
}

__global__ void SweY(double *Hs_in, double *h, double *Hs_out, double *K2n,
    double *K2s, double *Ztopo, double *mann, double *eRF, double*IN, 
    double *ET, int M, int N, int t) {
  int i0 = threadIdx.x + blockIdx.x * blockDim.x;
  int i = i0+1;

  // Declaration of local variable in each thread
  double a_y[M0-2], b_y[M0-2], c_y[M0-2], d_y[M0-2], f_y[M0-2];
  double c_star_y[M0-2], d_star_y[M0-2];

  while (i > 0 && i < N-1) {
    // Set the tri-diagonal systems in Y-direction
    SweYMatrix(a_y, b_y, c_y, d_y, M, N, i, t, Hs_in, h, mann, eRF, IN, ET, K2n, K2s);

    // Solve the tri-diagonal systems
    ThomasAlgorithm(a_y, b_y, c_y, d_y, f_y, c_star_y, d_star_y, M-2);

    // Update global values from TDMA solution
    for (int j = 1; j < M-1; j++) {
      Hs_out[j*N+i] = f_y[j-1];
      h[j*N+i] = Hs_out[j*N+i] - Ztopo[j*N+i];
    }
    i += gridDim.x * blockDim.x;
  }
}


// --------------------------------------------------------------------
// SWE_Explicit()
//    Solve the SWE for 2D overland flow using EXPLICIT method.
//    Update is implemeted via pointer
//
//    References:
//        * Wasantha Lal 1998, J. Hydraulic Engineering
//        * Shen & Phanikumar 2010, Adv. Water Resour.
//
//    Input...
//        * Hs_in: Water elevation [L]
//        * Ztopo: Ground surface elevation [L]
//        * h: Water depth [L]
//        * mann: Manning's coefficient [L^{-1/3} T]
//        * eRF: effective rainfall [L/T]
//        * ET: Evaporation [L/T]
//        * Scalar indexes
//
//    Output...
//        Matrix A and vector RHS d
//        * Hs_out: Water elevation[L]
//        * h: Water depth [L]
//        * IN: Infiltration [L/T]
//        * K2n and K2s: K values in north and south directions
//        * K2w and K2e: K values in east and west directions
// --------------------------------------------------------------------
__global__ void SWE_Explicit(double *Hs_in, double *h, double *Hs_out,
                  double *K2w, double *K2e, double *K2n, double *K2s,
                  double *Ztopo, double *mann, double *eRF, double*IN,
                  double *ET, int M, int N, int t) {
  int i0 = blockIdx.x * blockDim.x + threadIdx.x;
  int j0 = blockIdx.y * blockDim.y + threadIdx.y;
  int i = i0+1;
  int j = j0+1;

  double Keast, Kwest, Knorth, Ksouth, Sse, Ssw, Ssn, Sss;
  double h_, h_p1, h_m1, h_pN, h_mN, mann_, mann_p1, mann_m1, mann_pN, mann_mN;
  double Hs_, Hs_p1, Hs_pN, Hs_pNp1, Hs_pNm1, Hs_m1, Hs_mN, Hs_mNp1, Hs_mNm1;
  int I;

  while (i > 0 && i < N-1 && j > 0 && j < M-1) {
    I = j*N+i;

    // Transfer data from global to local memory for speed
    // p1: +1, m1: -1, pN: +N, mN: -N
    
    // water depth h
    h_ = h[I];
    h_p1 = h[I+1];
    h_pN = h[I+N];
    h_m1 = h[I-1];
    h_mN = h[I-N];

    // water elevation Hs=h+z
    Hs_ = Hs_in[I];
    Hs_p1 = Hs_in[I+1];
    Hs_pN = Hs_in[I+N];
    Hs_pNp1 = Hs_in[I+N+1];
    Hs_pNm1 = Hs_in[I+N-1];
    Hs_m1 = Hs_in[I-1];
    Hs_mN = Hs_in[I-N];
    Hs_mNp1 = Hs_in[I-N+1];
    Hs_mNm1 = Hs_in[I-N-1];

    // Manning's coefficients
    mann_ = mann[I];
    mann_p1 = mann[I+1];
    mann_pN = mann[I+N];
    mann_m1 = mann[I-1];
    mann_mN = mann[I-N];

    // Calculate slope Ss east and K east - Eqn (12), Lal 1998
    Sse = sqrt(pow((Hs_p1 - Hs_)/dx, 2.0) 
        + pow((Hs_pNp1+Hs_pN-Hs_mNp1-Hs_mN)/(4*dy), 2.0));

    if (Hs_p1 > Hs_) {
      if (h_p1 > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1+h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }
    } else {
      if (h_ > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1 + h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }
    }
    K2e[I] = Keast;

    // Calculate slope Ss west and K west - Eqn (12), Lal 1998
    Ssw = sqrt(pow((Hs_ - Hs_m1)/dx, 2.0)
        + pow((Hs_pN+Hs_pNm1-Hs_mN-Hs_mNm1)/(4*dy), 2.0));

    if (Hs_m1 > Hs_) {
      if (h_m1 > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }
    } else {
      if (h_ > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }
    }
    K2w[I] = Kwest;

    // Calculate slope Ss north and K north - Eqn (12), Lal 1998
    Ssn = sqrt(pow((Hs_pN - Hs_)/dy, 2.0)
        + pow((Hs_pNp1+Hs_p1-Hs_pNm1-Hs_m1)/(4*dx), 2.0));

    if (Hs_pN > Hs_) {
      if (h_pN > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }
    } else {
      if (h_ > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }
    }
    K2n[I] = Knorth;

    // Calculate slope Ss south and K south - Eqn (12), Lal 1998
    Sss = sqrt(pow((Hs_ - Hs_mN)/dy, 2.0)
        + pow((Hs_p1+Hs_mNp1-Hs_m1-Hs_mNm1)/(4*dx), 2.0));

    if (Hs_mN > Hs_) {
      if (h_mN > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }
    } else {
      if (h_ > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }
    }
    K2s[I] = Ksouth;

    // Update water elevation values via EXPLICIT scheme.
    Hs_out[I] = Hs_ + dt * (Kwest*(Hs_m1-Hs_)+Keast*(Hs_p1-Hs_)
        + Knorth*(Hs_pN-Hs_)+Ksouth*(Hs_mN-Hs_))/(dx*dy) + (eRF[I]-IN[I]-ET[t]);
    h[I] = Hs_out[I] - Ztopo[I];

    __syncthreads();
    i += gridDim.x * blockDim.x;

    if (i >= N-1) {
      j += gridDim.y * blockDim.y;
      i = i0+1;
    }
  }
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------