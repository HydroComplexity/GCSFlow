#include "../include/class.h"
#include "../include/variableclass.h"

void SimulationPerStep(ProjectClass *Project, OverlandClass *Swe2D,
		HostVarTime *HVTime,int M, int N, int P, int tt, int BSZ, int TSZ,
		dim3 dimGrid, dim3 dimBlockX, dim3 dimBlockY, dim3 dimBlockZ);

void SaveOutlet(FileNameClass *File, OverlandClass *Swe2D, HostVarTime *HVTime, 
		int M, int N, int P, int time_data, int time_steps);

void SavePerStep(ProjectClass *Project, HostData2D *Host2D, HostData3D *Host3D,
    FileNameClass *File, int M, int N, int P, int t_data);
