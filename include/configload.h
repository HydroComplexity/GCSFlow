#ifndef LOADCONFIGFILE_H
#define LOADCONFIGFILE_H

#include "../include/class.h"
void LoadConfigFile(const char *file_config, ProjectClass *Project,
    FileNameClass *File, DomainClass *Topo, SoilPropertiesClass *Soil,
    SubsurfaceClass *Richard3D, OverlandClass *Swe2D);
#endif      