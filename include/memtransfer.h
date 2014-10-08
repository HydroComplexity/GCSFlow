#include "../include/class.h"
#include "../include/variableclass.h"  

void cudaConstTransfer(ProjectClass *Project, DomainClass *Topo,
    SoilPropertiesClass *Soil, SubsurfaceClass *Richard3D, OverlandClass *Swe2D);

void TransferData(int M, int N, int P, int time_data, int time_steps,
    int substeps, HostData2D *Host2D, HostData3D *Host3D, HostDataTime *HDTime, 
    HostVarTime *HVTime);
