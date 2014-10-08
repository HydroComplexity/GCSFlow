#include "../include/class.h"

#ifndef NUMERICALMODEL_H
#define NUMERICALMODEL_H
int NumericalModel(ProjectClass *Project, FileNameClass *File,
    DomainClass *Topo, SoilPropertiesClass *Soil, SubsurfaceClass *Richard3D, 
    OverlandClass *Swe2D, HostData2D *Host2D, HostData3D *Host3D, 
    HostDataTime *HDTime, HostVarTime *HVTime);
#endif
