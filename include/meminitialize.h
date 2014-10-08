#include "../include/class.h"
#include "../include/variableclass.h"  

void cudaInitializeData(int M, int N, int P, int time_data, int time_steps);

void InitializeVariables(int M, int N, int P, int BSZ, int TSZ);

void FreeGPUMemory();