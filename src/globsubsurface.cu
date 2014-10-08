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


// This file provides global functions for subsurface flow model. Global 
// functions are linked with device functions. Both shared and non_shared 
// memory approaches are implemented.

#include "../include/constants.h"
#include "../include/devmaxfind.h"
#include "../include/devtdma.h"
#include "../include/devvanGenuchten.h"
#include "../include/devsubsurface.h"
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// vanGenuchten_Intial()
//    Estimates theta, K, and C from h using closed form formulation by
//    van Genuchten (1980).
//    Input: Vector pressure head h.
//    Output: Vector theta, K.
// --------------------------------------------------------------------
__global__ void vanGenuchtenIntial(double *theta, double *K, double *Ksat,
    double *h, int size) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  while (i < size) {
    double Se, _theta, _h, n, m, lambda;
    n = nv;
    lambda = n - 1.0;
    m = lambda/n;

    // Convert unit from [m] to [cm]
    _h = h[i] * 100;

    // Compute the volumetric moisture content [eqn 21]
    if (_h < 0) {
      _theta = (theta_S - theta_R) / pow(1.0 + pow((alpha*(-_h)), n), m) + theta_R;
    } else {
      _theta = theta_S;
    }
    theta[i] = _theta;

    // Compute the effective saturation [eqn 2]
    Se = (_theta - theta_R)/(theta_S - theta_R);

    // Compute the hydraulic conductivity [eqn 8] - Convert to unit: m/hr
    K[i] = Ksat[i] * sqrt(Se) * (1.0 - pow(1.0-pow(Se, 1.0/m), m))
        *(1.0 - pow(1.0 - pow(Se, 1.0/m), m));
           
    // Update threads if vector is long
    i += blockDim.x * gridDim.x;
  }
}


// --------------------------------------------------------------------
// RichZ()
// RichX()
// RichY()
//    These functions sets up the discritized equation for solving tridiagonal
//    systems in 3D subsurface flow model using ADI method.
//    RichZ() for vertical z-direction modeling.
//    RichX() and RichY() for horizontal x and y direction modeling.
//    References:
//        * Celia et al (1990) WRR,
//        * Clement et al (1994) J. Hydrol.,
//        * Sulis et al (2010) AWR.
//
//    Input...
//        * h_in: input pressure head 3D matrix [L]
//        * theta_in: input soil moisture 3D matrix [-]
//        * K_in: input hydraulic conductivity 3D matrix [L/T]
//        * h2D: input water depth 2D matrix [L]
//        * PPT: Rainfall [L/T]
//
//    Output...
//        Matrix A and vector RHS d
//        * h_out: out pressure head 3D matrix [L]
//        * theta_out: output soil moisture 3D matrix [-]
//        * K_out: output hydraulic conductivity 3D matrix [L/T]
// --------------------------------------------------------------------
__global__ void RichZ(double *h_in, double *theta_in, double *K_in,
    double *h_out, double *theta_out, double *K_out, double *Ksat, double *h2D,
    double PPT, double *IN, double *Psidiff, int *iter_z, int M, int N, int P) {
  // 2D thread and block id
  int i0 = blockIdx.x * blockDim.x + threadIdx.x;
  int j0 = blockIdx.y * blockDim.y + threadIdx.y;
  int i = i0+1;
  int j = j0+1;

  // Initialize local variables in z-direction
  double a_z[P0], b_z[P0], c_z[P0], d_z[P0], f_z[P0];
  double c_star_z[P0], d_star_z[P0];
  double hnp1m_z[P0], hnp1mp1_z[P0], Cnp1m_z[P0], Knp1m_z[P0], Ksat_z[P0];
  double thetan_z[P0], thetanp1m_z[P0], thetanp1mp1_z[P0];
  double htop, ph, qn;
  int TopBound, I;

  // Loop until thread id is out of the domain
  while (i > 0 && i < N-1 && j > 0 && j < M-1) {
    for (int k = 0; k < P; k++) {
      I = k*M*N + j*N + i;
      hnp1m_z[k] = h_in[I];
      thetan_z[k] = theta_in[I];
      Ksat_z[k] = Ksat[I];
    }

    // Determine top boundary
    ph = h2D[j*N + i];
    htop = ph;
    qn = PPT;

    if (ph > 0) {
      TopBound = 0;   // Dirichlet (head) boundary
      htop = ph;
    } else {
      TopBound = 1;   // Neumann (flux) boundary
    }

    int stop_flag = 0;
    int niter = 0;

    // Loop for convergence test
    while (stop_flag == 0 && niter < maxiterZ) {
      // Get C, K, theta from vanGenuchten function
      vanGenuchten(Cnp1m_z, Knp1m_z, Ksat_z, thetanp1m_z, hnp1m_z, nv, P);

      // Set up the tri-digonal system in z-direction
      RichZMatrix(a_z, b_z, c_z, d_z, Cnp1m_z, Knp1m_z, K_in, hnp1m_z, h_in,
          thetanp1m_z, thetan_z, i, j, M, N, P);

      // Solve the tri-diagonal systems using TDMA
      ThomasAlgorithm(a_z, b_z, c_z, d_z, f_z, c_star_z, d_star_z, P);
      niter += 1;

      // If error < stop tolerance, stop and move to next step
      if (MaxErrorZ(f_z, P) < stop_tol) {
        stop_flag = 1;
        for (int k = 0; k < P; k++) {
          hnp1mp1_z[k] = hnp1m_z[k] + f_z[k];
        }

        // Force top boundary conditions
        if (TopBound == 0) {
          hnp1mp1_z[0] = htop;
        } else {
          hnp1mp1_z[0] = hnp1mp1_z[1] - dz + 2*dz*qn/dt/(Knp1m_z[0]+Knp1m_z[1]);
        }

        // Force bottom boundary conditions
        hnp1mp1_z[P-1] = hnp1mp1_z[P-2] + dz;

        // Get C, K, theta from vanGenuchten function
        vanGenuchten(Cnp1m_z, Knp1m_z, Ksat_z, thetanp1m_z, hnp1mp1_z, nv, P);

        // Update solution from vectors to the global matrix
        for (int k = 0; k < P; k++) {
          theta_out[k*M*N+j*N+i] = thetanp1m_z[k];
          thetanp1mp1_z[k] = thetanp1m_z[k];
          K_out[k*M*N+j*N+i] = Knp1m_z[k];
          h_out[k*M*N+j*N+i] = hnp1mp1_z[k];
        }
      } else {    // If error > stop tolerance, run the loop again.
        for (int k = 0; k < P; k++) {
          if (abs(f_z[k]) > 2.0)
            f_z[k] = 0.001;
          hnp1mp1_z[k] = hnp1m_z[k] + f_z[k];
          hnp1m_z[k] = hnp1mp1_z[k];
        }

        // Force boundary conditions
        if (TopBound == 0) {
          hnp1m_z[0] = htop;
        } else {
          hnp1m_z[0] = hnp1m_z[1] - dz + 2*dz*qn/dt/(Knp1m_z[0]+Knp1m_z[1]);
        }

        // Boundary switching based on top node's state
        // Switch to Dirichlet boundary
        if (hnp1m_z[0] >= psi_min) {
          TopBound = 0;
          if (ph > 0) {
            htop = ph;
          } else {
            htop = psi_min;
          }

        } else {    // Switch to Neumann boundary
          TopBound = 1;
          qn = PPT;
        }

        // Force bottom boundary conditions
        if (BoundBottom == 0){
          hnp1m_z[P-1] = Psi_Bottom;
        } else {
          hnp1m_z[P-1] = hnp1m_z[P-2] + dz;
        }
      }
    }

    iter_z[j*N+i] = niter;

    // Update pressure head, moisture and calculate Infiltration
    if (niter < maxiterZ) {
      for (int k = 0; k < P; k++) {
        h_out[k*M*N+j*N+i] = hnp1mp1_z[k];
        theta_out[k*M*N+j*N+i] = thetanp1mp1_z[k];
      }
      IN[j*N+i] = -(Knp1m_z[0]) * (hnp1mp1_z[1] - hnp1mp1_z[0] - dz)/dz*dt;
      //Psidiff[j*N+i] = hnp1mp1_z[1] - hnp1mp1_z[0]; 
      Psidiff[j*N+i] = (hnp1mp1_z[0]);
    } else {
      for (int k = 0; k < P; k++) {
        h_out[k*M*N+j*N+i] = hnp1m_z[k];
        theta_out[k*M*N+j*N+i] = thetanp1m_z[k];
      }
      IN[j*N+i] = -(Knp1m_z[0]) * (hnp1m_z[1] - hnp1m_z[0] - dz)/dz*dt;
      //Psidiff[j*N+i] = hnp1m_z[1] - hnp1m_z[0];
      Psidiff[j*N+i] = (hnp1mp1_z[0]);
    }

    __syncthreads();    // All thread must sync at this point
    i += gridDim.x * blockDim.x;

    if (i >= N-1) {
      j += gridDim.y * blockDim.y;
      i = i0+1;
    }
  }
}

__global__ void RichX(double *h_in, double *theta_in, double *K_in,
    double *h_out, double *theta_out, double *K_out, double *Ksat, int * iter_x,
    int M, int N, int P) {
  // 2D thread and block id
  int j0 = blockIdx.x * blockDim.x + threadIdx.x;
  int k0 = blockIdx.y * blockDim.y + threadIdx.y;
  int j = j0+1;
  int k = k0+1;

  // Initialize local variables in x-direction
  double a_x[N0-2], b_x[N0-2], c_x[N0-2], d_x[N0-2], f_x[N0-2];
  double c_star_x[N0-2], d_star_x[N0-2];
  double hnp1m_x[N0], hnp1mp1_x[N0], Cnp1m_x[N0], Knp1m_x[N0], Ksat_x[N0];
  double thetan_x[N0], thetanp1m_x[N0];
  int I;

  while (j > 0 && j < M-1 && k > 0 && k < P-1) {
    for (int i = 0; i < N; i++) {
      I = k*M*N + j*N + i;
      hnp1m_x[i] = h_in[I];
      thetan_x[i] = theta_in[I];
      Ksat_x[i] = Ksat[I];      
    }

    int stop_flag = 0;              // Define a dummy stopping variable
    int niter = 0;                  // Define an iteration counter

    // Loop for convergence test
    while (stop_flag == 0 && niter < maxiterX) {
      // Get C, K, theta from vanGenuchten function
      vanGenuchten(Cnp1m_x, Knp1m_x, Ksat_x, thetanp1m_x, hnp1m_x, nv, N);

      // Set up the tri-digonal system in x-direction
      RichXMatrix(a_x, b_x, c_x, d_x, Cnp1m_x, Knp1m_x, K_in, hnp1m_x, h_in,
          thetanp1m_x, thetan_x, j, k, M, N, P);

      // Solve the tri-diagonal systems using TDMA
      ThomasAlgorithm(a_x, b_x, c_x, d_x, f_x, c_star_x, d_star_x, N-2);
      niter += 1;

      // If error < stop tolerance, stop and move to next step
      if (MaxError(f_x, N-2) < stop_tol) {
        stop_flag = 1;

        for (int i = 1; i < N-1; i++) {
          hnp1mp1_x[i] = hnp1m_x[i] + f_x[i-1];
        }

        // Force boundary conditions
        hnp1mp1_x[0] = hnp1mp1_x[1];
        hnp1mp1_x[N-1] = hnp1mp1_x[N-2];

        // Get C, K, theta from vanGenuchten function
        vanGenuchten(Cnp1m_x, Knp1m_x, Ksat_x, thetanp1m_x, hnp1mp1_x, nv, N);

        // Update solution from vectors to the global matrix
        for (int i = 0; i < N; i++) {
          theta_out[k*M*N+j*N+i] = thetanp1m_x[i];
          K_out[k*M*N+j*N+i] = Knp1m_x[i];
          h_out[k*M*N+j*N+i] = hnp1mp1_x[i];
        }

      } else {    // If error > stop tolerance, run the loop again.
        for (int i = 1; i < N-1; i++) {
            hnp1mp1_x[i] = hnp1m_x[i] + f_x[i-1];
            hnp1m_x[i] = hnp1mp1_x[i];
        }

        // Force boundary conditions
        if (BoundEast == 0){
          hnp1m_x[0] = Psi_East;
        } else {
          hnp1m_x[0] = hnp1m_x[1];
        }
        
        if (BoundWest == 0){
          hnp1m_x[N-1] = Psi_West;  
        } else {
          hnp1m_x[N-1] = hnp1m_x[N-2];
        }
      }
    }
    iter_x[k*M+j] = niter;        
    __syncthreads();

    j += gridDim.x * blockDim.x;
    if (j >= M-1) {
      k += gridDim.y * blockDim.y;
      j = j0+1;
    }
  }
}

__global__ void RichY(double *h_in, double *theta_in, double *K_in,
    double *h_out, double *theta_out, double *K_out, double *Ksat, int *iter_y,
    int M, int N, int P) {
  // 2D thread and block id
  int i0 = blockIdx.x * blockDim.x + threadIdx.x;
  int k0 = blockIdx.y * blockDim.y + threadIdx.y;
  int i = i0+1;
  int k = k0+1;

  // Initialize local variables in y-direction
  double a_y[M0-2], b_y[M0-2], c_y[M0-2], d_y[M0-2], f_y[M0-2];
  double c_star_y[M0-2], d_star_y[M0-2];
  double hnp1m_y[M0], hnp1mp1_y[M0],  Cnp1m_y[M0], Knp1m_y[M0], Ksat_y[M0];
  double thetan_y[M0], thetanp1m_y[M0];
  int I;

  while (i > 0 && i < N-1 && k > 0 && k < P-1) {
    for (int j = 0; j < M; j++) {
      I = k*M*N + j*N + i;
      hnp1m_y[j] = h_in[I];
      thetan_y[j] = theta_in[I];
      Ksat_y[j] = Ksat[I];
    }
    int stop_flag = 0;              // Define a dummy stopping variable
    int niter = 0;                  // Define an iteration counter

    // Loop for convergence test
    while (stop_flag == 0 && niter < maxiterY) {
      // Get C, K, theta from vanGenuchten function
      vanGenuchten(Cnp1m_y, Knp1m_y, Ksat_y, thetanp1m_y, hnp1m_y, nv, M);

      // Set up the tri-digonal system in y-direction
      RichYMatrix(a_y, b_y, c_y, d_y, Cnp1m_y, Knp1m_y, K_in, hnp1m_y, h_in,
          thetanp1m_y, thetan_y, i, k, M, N, P);

      // Solve the tri-diagonal systems using TDMA
      ThomasAlgorithm(a_y, b_y, c_y, d_y, f_y, c_star_y, d_star_y, M-2);

      niter += 1;

      if (MaxError(f_y, M-2) < stop_tol) {
        stop_flag = 1;
        for (int j = 1; j < M-1; j++) {
          hnp1mp1_y[j] = hnp1m_y[j] + f_y[j-1];
        }

        // Force boundary conditions
        hnp1mp1_y[0] = hnp1mp1_y[1];
        hnp1mp1_y[M-1] = hnp1mp1_y[M-2];

        // Get C, K, theta from vanGenuchten function
        vanGenuchten(Cnp1m_y, Knp1m_y, Ksat_y, thetanp1m_y, hnp1mp1_y, nv, M);

        // Update solution from vectors to the global matrix
        for (int j = 0; j < M; j++) {
          theta_out[k*M*N+j*N+i] = thetanp1m_y[j];
          K_out[k*M*N+j*N+i] = Knp1m_y[j];
          h_out[k*M*N+j*N+i] = hnp1mp1_y[j];
        }
      } else {
        for (int j = 1; j < M-1; j++) {
          hnp1mp1_y[j] = hnp1m_y[j] + f_y[j-1];
          hnp1m_y[j] = hnp1mp1_y[j];
        }

        // Force boundary conditions
        if (BoundNorth == 0){
          hnp1m_y[0] = Psi_North;
        } else {
          hnp1m_y[0] = hnp1m_y[1];
        }

        if (BoundSouth == 0){
          hnp1m_y[M-1] = Psi_South;
        } else {
          hnp1m_y[M-1] = hnp1m_y[M-2];
        }
      }
    }
    iter_y[k*N+i] = niter;            
    __syncthreads();

    i += gridDim.x * blockDim.x;
    if (i >= N-1) {
      k += gridDim.y * blockDim.y;
      i = i0+1;
    }
  }
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------
