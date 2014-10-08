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


// This file provide a function to print information loaded from config file
// Special classes defined in class.h and functions in parser.h are used.

#include <stdio.h>
#include "../include/class.h"
#include "../include/colorcode.h"

// --------------------------------------------------------------------
// PrintModelParameters()
//    Print out Project information and Model parameters
//
//    - dt_h: time step in host memory
//    - time_steps: Number of timesteps
// --------------------------------------------------------------------
void PrintModelInfo(ProjectClass *Project, SubsurfaceClass *Richard3D, 
    DomainClass *Topo, OverlandClass *Swe2D, SoilPropertiesClass *Soil) {
  
  double dt_h = Project->dt_data/Project->substeps;
  int time_steps = Project->time_steps;

  printf("\n");
  printf(L1 " GPU-BASED CONJUNCTIVE SURFACE-SUBSURFACE FLOW MODEL \n" RS);
  
  printf(L2"   1. Name: " RS "%s \n", Project->Name);

  printf(L2"   2. Domain Size: \n" RS);
  printf("\t # of cells in x-direction: %d \n", Topo->Nx);
  printf("\t # of cells in y-direction: %d \n", Topo->My);
  printf("\t # of cells in z-direction: %d \n", Topo->Pz);
  printf("\n");

  printf(L2"   3. Spatial Resolutions: \n" RS);
  printf(" \t dx = %5.2f [m]\n", Topo->GridSizeDx);
  printf(" \t dy = %5.2f [m]\n", Topo->GridSizeDy);
  printf(" \t dz = %5.2f [m]\n", Topo->GridSizeDz);
  printf("\n");

  printf(L2"   4. Simulation Time: " RS "%4.2f [hr] \n", dt_h*time_steps);
  printf("\n");

  printf(L2"   5. Time steps:\n" RS);
  printf("\t Data: %-6.4f [hr] \n", Project->dt_data);
  printf("\t Computation: %-6.4f [hr] \n" RS, dt_h);
 	printf("\n");

  printf(L1 " MODEL PARAMETERS:\n" RS);
	printf(L2"   1. Overland Flow \n" RS);
  printf("\t Delta = %10.8f\n", Swe2D->Delta);
	printf("\t Minimum Depth = %6.4f [m]\n", Swe2D->MinimumDepth);
	printf("\t Critical Depth = %4.2f [m]\n", Swe2D->CriticalDepth);
	printf("\t K0 = %4.2f\n", Swe2D->K0);
	printf("\n");

	printf(L2"   2. Subsurface Flow \n" RS);
  printf("\t Minimum Pressure Head = %5.3f [m]\n", Richard3D->Psimin);
	printf("\t Stopping Tolerance = %6.4f\n", Richard3D->stoptol);
	printf("\t Alpha = %3.2f\n", Soil->Alpha);
	printf("\t Pore-size distribution = %3.2f\n", Soil->n);
	printf("\t Porosity = %3.2f\n", Soil->porosity);
	printf("\t Residual water content = %4.2f\n", Soil->Theta_R);
	printf("\t Saturated water content = %4.2f\n", Soil->Theta_S);
  printf("\n"); printf("\n"); 

  printf(L1 " SIMULATION PROGRESS...\n" RS);
}

