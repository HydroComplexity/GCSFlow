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


// This file declares __constant__ and pointer variables used in device.

// Grid cells and time resolution
__constant__ double dx;             // horizontal x direction [m]
__constant__ double dy;             // horizontal y direction [m]
__constant__ double dz;             // vertical z direction [m]
__constant__ double dt;             // time step [hour]

// Surface model parameters
__constant__ double delta;          // [-]
__constant__ double hmin;           // [m]
__constant__ double hcri;           // Critical depth [m]
__constant__ double K0;             // [m2/h]

// Subsurface model parameters
__constant__ double alpha;          // [1/cm]
__constant__ double theta_S;        // [cm3/cm3]
__constant__ double theta_R;        // [cm3/cm3]
__constant__ double lambda;         // [-]
__constant__ double Ss;             // [1/m]
__constant__ double poros;          // [-]
__constant__ double nv;             // [-]

__constant__ double psi_min;        // [m]
__constant__ double stop_tol;       // [m]
__constant__ int BoundBottom;       // [-]
__constant__ int BoundNorth;        // [-]
__constant__ int BoundSouth;        // [-]
__constant__ int BoundEast;         // [-]
__constant__ int BoundWest;         // [-]
__constant__ int maxiterZ;          // [-]
__constant__ int maxiterX;          // [-]
__constant__ int maxiterY;          // [-]

__constant__ double Psi_Bottom;     // [m]
__constant__ double Psi_North;      // [m]
__constant__ double Psi_South;      // [m]
__constant__ double Psi_East;       // [m]
__constant__ double Psi_West;       // [m]

//-----------------------------------------------------------------------
// Variable declaration
// type *name           // Description - [dimensions] - [unit]
//-----------------------------------------------------------------------
double *Psi_in_d;       // pressure head input - [M,N,P] - [m]
double *Psi_out_d;      // pressure head output - [M,N,P] - [m]
double *theta_in_d;     // soil moisture input - [M,N,P] - [m]
double *theta_out_d;    // soil moisture output - [M,N,P] - [m]
double *K_in_d;         // hydraulic conductivity input - [M,N,P] - [m]
double *K_out_d;        // hydraulic conductivity output - [M,N,P] - [m]
double *Ksat_d;         // Saturated hydraulic conductivity - [M,N,P] - [m]
double *C_in_d;         // specific moisture storage input - [M,N,P] - [1/m]
double *Hs_in_d;        // water elevation input - [M,N] - [m]
double *Hs_out_d;       // water elevation output - [M,N] - [m]
double *h_d;            // water depth - [M,N] - [m]
double *eRF;            // effective rainfall - [M,N] - [m]
double *IN;             // infiltration - [M,N] - [m]
double *Psidiff;        // Head differece between 1st and 2nd layers - [M,N]
int *iter_z_d;
int *iter_x_d;
int *iter_y_d;
double *htop;           //
double *hbottom;        //
double *PPT_d;          // precipitation - [Time] -[m]
double *ET_d;           // evaporation - [Time] -[m]
double *Ztopo_d;        // topography - [M,N] - [m]
double *mann_d;         // Manning's roughness coefficient - [M,N] - [m]
double *K2w_d;          // K2D in west - [M,N] - [m]
double *K2e_d;          // K2D in east - [M,N] - [m]
double *K2n_d;          // K2D in north - [M,N] - [m]
double *K2s_d;          // K2D in south - [M,N] - [m]
double *u2_d;           // velocity in x direction - [M,N] - [m]
double *v2_d;           // velocity in y direction - [M,N] - [m]
double *hout_d;         // water depth at outlet - [Time] -[m]
double *vout_d;         // u velocity at outlet - [Time] -[m]
double *uout_d;         // v velocity at outlet - [Time] -[m]
double *qss_d;          // discharge at outlet - [Time] -[m]
