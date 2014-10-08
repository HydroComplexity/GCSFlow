#include "../include/class.h"
#include "../include/variableclass.h"

#ifndef GETDATAINFO_H
#define GETDATAINFO_H
  void GetDataInfo(ProjectClass *Project, FileNameClass *File,
      DomainClass *Topo, DimensionClass *Dimension);
#endif
                                    
#ifndef GETFILE_H
#define GETFILE_H
  void GetFile(const char *file_name, const char *var_name, int ndims, int *dim);
#endif

#ifndef LOADDATATOHOST_H
#define LOADDATATOHOST_H
  void LoadDataToHost(FileNameClass *File, DimensionClass *Dimension, 
      HostData2D *Host2D, HostData3D *Host3D, HostDataTime *HDTime);
#endif

#ifndef LOADFILE_H
#define LOADFILE_H
  void LoadFile(const char *file_name, const char *var_name, int ndims, 
      int *dim, double *data);
#endif

