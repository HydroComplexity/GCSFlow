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


// In this file, we provide functions to check and load input data files.
// A NetCDF file is first examine to identify the dimensions of the data
// using GetDataInfo(). The length of data inside NetCDF files are
// evaluated using GetFile(). Second, The data is loaded into the host
// memory using LoadDataToHost() and LoadFile() functions.
// Special classes defined in class.h and varclass.h are used.

#include <netcdf.h>
#include "../include/class.h"
#include "../include/variableclass.h"


// --------------------------------------------------------------------
// GetDataInfo()
//    Gets info about a NetCDF input file. This returns the number of
//    dimension of a variable and length of each dimension.
//
//    - ncid: NetCDF ID
//    - var_id: Variable ID
//    - dimids: Dimension IDs
//    - file_name: NetCDF input file
//    - dim: dimensions of the dataset
// --------------------------------------------------------------------
void GetFile(const char *file_name, const char *var_name, int ndims, int *dim) {
  int ncid, var_id, dimids[3];
  size_t length;

  // Open Netcdf file with NC_NOWRITE options
  nc_open(file_name, NC_NOWRITE, &ncid);

  // Get variable id, dimension, size
  nc_inq_varid(ncid, var_name, &var_id);
  nc_inq_vardimid(ncid, var_id, dimids);

  for (int i = 0; i < ndims; i++) {
    nc_inq_dimlen(ncid, dimids[i], &length);
    dim[i] = length;
  }

  // Close Netcdf file
  nc_close(ncid);
}


// --------------------------------------------------------------------
// GetDataInfo()
//    Gets info about a data input files. Use GetFile() function to
//    get Dimensions and size of data files.
//    
//    - Ztopo: Topographic elevation [L]
//    - mann: Manning's coefficient
//    - Psi_init: Initial pressure head values
//    - PPT_in: Input precipitation from data file
//    - Nx, My, Pz: Number of cells in x, y, and z directions, respectively
//    - dim_XXXX: dimension of variable XXXX
// --------------------------------------------------------------------
void GetDataInfo(ProjectClass *Project, FileNameClass *File,
    DomainClass *Topo, DimensionClass *Dimension) {
  GetFile(File->topo_file, "Ztopo", 2, Dimension->dim_topo);
  GetFile(File->parameter_file, "mann", 2, Dimension->dim_para);
  GetFile(File->parameter_file, "Ksat", 3, Dimension->dim_Ksat);  
  GetFile(File->init_conds_file, "Psi_init", 3, Dimension->dim_Psi);
  GetFile(File->forcing_file, "PPT_in", 1, Dimension->dim_forcing);

  // Get and check domain dimensions
  Topo->Pz = Dimension->dim_Psi[0];
  Topo->My = Dimension->dim_Psi[1];
  Topo->Nx = Dimension->dim_Psi[2];

  // Get number of time steps for simulation
  Project->time_data = Dimension->dim_forcing[0];
  Project->time_steps = Project->time_data * Project->substeps;
}


// --------------------------------------------------------------------
// LoadData()
//    Loads data from a variable of a NetCDF input file. This transfers
//    the data via pointer. Works only for double type.
//
//    - ncid: NetCDF ID
//    - var_id: Variable ID
//    - *file_name: NetCDF input file
//    - *var_name: variable name in NetCDF file
//    - *dim: dimensions of the dataset
//    - *data: output data loaded from NetCDF file to CPU memory
// --------------------------------------------------------------------
void LoadFile(const char *file_name, const char *var_name, int ndims,
    int *dim, double *data) {
  int ncid, var_id;

  // Open Netcdf file with NC_NOWRITE options
  nc_open(file_name, NC_NOWRITE, &ncid);

  // Get variable id
  nc_inq_varid(ncid, var_name, &var_id);

  // Get data based on data_type
  nc_get_var_double(ncid, var_id, &data[0]);

  // Close Netcdf file
  nc_close(ncid);
}


// --------------------------------------------------------------------
// LoadDataToHost()
//    Load input data from NetCDF files to Host memory.
// --------------------------------------------------------------------
void LoadDataToHost(FileNameClass *File, DimensionClass *Dimension,
          HostData2D *Host2D, HostData3D *Host3D, HostDataTime *HDTime) {
  // Load data from NetCDF files
  LoadFile(File->topo_file, "Ztopo", 2, Dimension->dim_topo, Host2D->Ztopo);
  LoadFile(File->init_conds_file, "Psi_init", 3, Dimension->dim_Psi, 
      Host3D->Psi_in);
  LoadFile(File->parameter_file, "Ksat_in", 3, Dimension->dim_Ksat, Host3D->Ksat);
  LoadFile(File->parameter_file, "mann", 2, Dimension->dim_para, Host2D->mann);  
  LoadFile(File->forcing_file, "PPT_in", 1, Dimension->dim_forcing, 
      HDTime->PPT_data);
  LoadFile(File->forcing_file, "Evap_in", 1, Dimension->dim_forcing, 
      HDTime->ET_data);
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------
