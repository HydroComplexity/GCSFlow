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


// In this file, we provide functions to load the configuration file.
// First, the configuration file from argv[1] is loaded into the program.
// The function scans everyline to load appropirate parameters/names
// for the model.
// Special classes defined in class.h and functions in parser.h are used.

#include "../include/class.h"
#include "../include/parser.h"

// --------------------------------------------------------------------
// LoadConfigFile()
//    Load info from configuration file input from argv[1]
//    Info is stored into classes defined in class header (class.h)
// --------------------------------------------------------------------
void LoadConfigFile(const char *file_config, ProjectClass *Project,
    FileNameClass *File, DomainClass *Topo, SoilPropertiesClass *Soil,
    SubsurfaceClass *Richard3D, OverlandClass *Swe2D) {
  
  // Parsing configuration file
  ParserConfigFile(file_config);

  // Load information/data into Project class
  Project->Name = GetOptionToChar("ProjectName");
  Project->BLOCK_DIMX = GetOptionToInt("BlockDimensionX");
  Project->BLOCK_DIMY = GetOptionToInt("BlockDimensionY");
  Project->BLOCK_DIMZ = GetOptionToInt("BlockDimensionZ");
  Project->BSZ = GetOptionToInt("BlockSize");
  Project->TSZ = GetOptionToInt("ThreadSize");  

  Project->PrintSubStep = GetOptionToInt("PrintPerSubStep");
  Project->SaveStep = GetOptionToInt("SavePerStep");
  Project->Save2DOverland = GetOptionToInt("Save2DOverland");
  Project->Save3DSubsurface = GetOptionToInt("Save3DSubsurface");
  Project->SaveIteration = GetOptionToInt("SaveIteration");
  Project->SaveInfiltration = GetOptionToInt("SaveInfiltration");
  Project->BC2D = GetOptionToInt("BoundaryCondition2D");  
  Project->dt_data = GetOptionToDouble("DataTimeStep");
  Project->substeps = GetOptionToInt("NumberSubTimeStep");
  
  // Load filename information into File class
  File->topo_file = GetOptionToChar("Topography");
  File->init_conds_file = GetOptionToChar("InitialCondition");
  File->forcing_file = GetOptionToChar("Forcing");
  File->parameter_file = GetOptionToChar("Parameters");
  File->output3D = GetOptionToChar("Output3D");
  File->output2D = GetOptionToChar("Output2D");
  File->output_veloc2D = GetOptionToChar("OutputVelocity2D");
  File->output_K2D = GetOptionToChar("OutputK2D");
  File->output_h = GetOptionToChar("OutputH");
  File->output_u = GetOptionToChar("OutputU");
  File->output_v = GetOptionToChar("OutputV");
  File->output_qss = GetOptionToChar("OutputQss");

  // Load Topographic information into Topo class
  Topo->GridSizeDx = GetOptionToDouble("GridSize.Dx");
  Topo->GridSizeDy = GetOptionToDouble("GridSize.Dy");
  Topo->GridSizeDz = GetOptionToDouble("GridSize.Dz");

  // Load vanGenuchten parameters into Soil class
  Soil->Alpha = GetOptionToDouble("Alpha");
  Soil->Theta_S = GetOptionToDouble("SaturatedWaterContent");
  Soil->Theta_R = GetOptionToDouble("ResidualWaterContent");
  Soil->n = GetOptionToDouble("PoreSizeDistribution");
  Soil->porosity = GetOptionToDouble("Porosity");
  Soil->Ss = GetOptionToDouble("SpecificStorage");
  
  // Load paramters to 3D subsurface flow model
  Richard3D->Psimin = GetOptionToDouble("MinimumPressureHead");
  Richard3D->stoptol = GetOptionToDouble("StopTolerance");
  Richard3D->BoundBottom = GetOptionToInt("BoundBottom");
  Richard3D->BoundNorth = GetOptionToInt("BoundNorth");
  Richard3D->BoundSouth = GetOptionToInt("BoundSouth");
  Richard3D->BoundEast = GetOptionToInt("BoundEast");
  Richard3D->BoundWest = GetOptionToInt("BoundWest");
  Richard3D->PsiBottom = GetOptionToDouble("PsiBottom");
  Richard3D->PsiNorth = GetOptionToDouble("PsiNorth");
  Richard3D->PsiSouth = GetOptionToDouble("PsiSouth");
  Richard3D->PsiWest = GetOptionToDouble("PsiWest");
  Richard3D->PsiEast = GetOptionToDouble("PsiEast");

  Richard3D->MaxIterX = GetOptionToInt("MaxIterationX");
  Richard3D->MaxIterY = GetOptionToInt("MaxIterationY");
  Richard3D->MaxIterZ = GetOptionToInt("MaxIterationZ");

  // Load parameters to 2D surface flow model parameters
  Swe2D->Numeric = GetOptionToDouble("NumericalScheme");
  Swe2D->Delta = GetOptionToDouble("Delta");
  Swe2D->MinimumDepth = GetOptionToDouble("MinimumDepth");
  Swe2D->CriticalDepth = GetOptionToDouble("CriticalDepth");
  Swe2D->K0 = GetOptionToDouble("K0");
}
