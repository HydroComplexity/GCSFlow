import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import os
from osgeo import gdal
import scipy.io as sio

# REPLACE YOUR TIF FILE FOR TOPOGRAPHY HERE
gdata = gdal.Open("Goose_Lidar/goose_rec.tif")
geo = gdata.GetGeoTransform()
data = gdata.ReadAsArray()

# Change the size of the study area here
Topography_input = data[0:100,0:200]

My,Nx = Topography_input.shape
Pz = 10
Ksat_input = 5e-3 * np.ones((My,Nx))
mann_input = 0.03 * np.ones((My,Nx))
Psi_init = -0.4 * np.ones((Pz,My,Nx))

# Check/create folder and topography file name
folder = 'output'

output_topo = 'topography_goosebw.nc'
output_para = 'parameters_goosebw.nc'
output_initconds = 'initial_conds_goosebw.nc'
if not os.path.exists(folder):
    os.makedirs(folder)



################################################################
# DO NOT CHANGE BELOW THIS LINE --------------------------------
################################################################

topography = Dataset(folder + '/' + output_topo, 'w', format='NETCDF4')
topography.description = 'Topography Data for GCSFlow'

para = Dataset(folder + '/' + output_para, 'w', format='NETCDF4')
para.description = 'Saturated hydraulic conductivity for GCSFlow'

initconds = Dataset(folder + '/' + output_initconds, 'w', format='NETCDF4')
initconds.description = 'Initial conditions for GCSFlow'

# Get dimension
My,Nx  = Topography_input.shape

# Set up dimension for files in NetCDF
topography.createDimension('y', My)
topography.createDimension('x', Nx)
para.createDimension('y', My)
para.createDimension('x', Nx)
initconds.createDimension('y', My)
initconds.createDimension('x', Nx)
initconds.createDimension('z', Pz)

"""
Create variables in the netcdf file
	var = netcdf.createVariable('Var_name', 'var_type', ('dimension_type'))
"""
Ztopo = topography.createVariable('Ztopo', 'f8', ('y','x'))
HK = para.createVariable('Ksat_in', 'f8', ('y','x'))
Mannings = para.createVariable('mann', 'f8', ('y','x'))
Psi_ic = initconds.createVariable('Psi_init', 'f8', ('z','y','x'))

# Assign data to variables in NetCDF file
Ztopo[:] = Topography_input
HK[:] = Ksat_input
Mannings[:] = mann_input
Psi_ic[:] = Psi_init
#HK[:] = 0.0025

# Close the file
topography.close()
para.close()
initconds.close()

# Plot topography for checking
fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))

im = ax0.imshow(Topography_input, cmap=plt.cm.jet)
fig.colorbar(im, ax=ax0)

im = ax1.imshow(Ksat_input, cmap=plt.cm.Spectral)
fig.colorbar(im, ax=ax1)

plt.show()