void SaveOutput2D(const char *file, int Nx, int My, const char *data_name1, 
    double *data_in1, const char *data_name2, double *data_in2);

void SaveOneIntOutput2D(const char *file, int Nx, int My, const char *data_name,
    int *data_in);
  
void SaveOneOutput2D(const char *file, int Nx, int My, const char *data_name, 
    double *data_in);

void SaveOutput3D(const char *file, int Nx, int My, int Pz, const char *data_name1, 
    double *data_in1, const char *data_name2, double *data_in2);

void SaveIterationMap(const char *file, int Nx, int My, int Pz, const char *nameZ, 
    int *Z_in, const char *nameX, int *X_in, const char *nameY, int *Y_in);

void SaveVariables1D(const char *file, int time_steps, const char *var_name, double *var);
  
void SaveK2D(const char *file, int My, int Nx, const char *nameW, double *data_W, 
    const char *nameE, double *data_E, const char *nameN, double *data_N, 
    const char *nameS, double *data_S);