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


#include <stdio.h>
#include "../include/gpuinfo.h"
#include "../include/dataload.h"
#include "../include/parser.h"
#include "../include/configload.h"
#include "../include/class.h"
#include "../include/variableclass.h"
#include "../include/modelkern.h"
#include "../include/printmodelinfo.h"

int main(int argc, char **argv) {
  int device[2];
  const char *file_config;

  // Variables with struct defined in class header
  ProjectClass *Project = new ProjectClass;
  OverlandClass *Swe2D = new OverlandClass;
  SubsurfaceClass *Richard3D = new SubsurfaceClass;
  SoilPropertiesClass *Soil = new SoilPropertiesClass;
  DomainClass *Topo = new DomainClass;
  FileNameClass *File = new FileNameClass;
  DimensionClass *Dimension = new DimensionClass;

  printf("\n STARTING... \n");

  // Load Configuration File
  if ( argc != 2 ) {
    printf("\t Error: No configuration input file!\n");
    printf("\t Usage: %s <configuration file>\n", argv[0]);
    printf("\t Exit!\n");
    exit(1);
  } else {
    file_config = argv[1];
  }
  LoadConfigFile(file_config, Project, File, Topo, Soil, Richard3D, Swe2D);

  // Get GPU device information & eligibility
  GPUGetInfo(device);

  // Get Information about input files (NetCDF)
  GetDataInfo(Project, File, Topo, Dimension);

  // Allocate 4 different types of variables in Host memory
  HostData2D *Host2D = new HostData2D(Topo->My*Topo->Nx, Topo->My*Topo->Pz, 
      Topo->Nx*Topo->Pz);
  HostData3D *Host3D = new HostData3D(Topo->My*Topo->Nx*Topo->Pz);
  HostDataTime *HDTime = new HostDataTime(Project->time_data);
  HostVarTime *HVTime = new HostVarTime(Project->time_steps);

  // Load Data from file to Host Memory
  LoadDataToHost(File, Dimension, Host2D, Host3D, HDTime);

  // Print out model parameters
  PrintModelInfo(Project, Richard3D, Topo, Swe2D, Soil);

  // Computation code runs in CUDA
  NumericalModel(Project, File, Topo, Soil, Richard3D, Swe2D, Host2D, Host3D,
      HDTime, HVTime);

  // Clear host memory
  delete[] Project; delete[] Swe2D; delete[] Richard3D; delete[] Soil; 
  delete[] Topo; delete[] File; delete[] Dimension; delete[] Host2D; 
  delete[] Host3D; delete[] HDTime; delete[] HVTime;
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------
