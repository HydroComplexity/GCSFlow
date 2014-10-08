__global__ void vanGenuchtenIntial(double *theta, double *K, double *Ksat, 
    double *h, int size);

__global__ void RichZ(double *h_in, double *theta_in, double *K_in, 
    double *h_out, double *theta_out, double *K_out, double *Ksat, double *h2D,
    double PPT, double *IN, double *Psidiff, int *iter_z_d, int M, int N, int P);

__global__ void RichX(double *h_in, double *theta_in, double *K_in, 
    double *h_out, double *theta_out, double *K_out, double *Ksat, 
    int *iter_x_d, int M, int N, int P);

__global__ void RichY(double *h_in, double *theta_in, double *K_in, 
    double *h_out, double *theta_out, double *K_out, double *Ksat, 
    int *iter_y_d, int M, int N, int P);