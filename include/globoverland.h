__global__ void SweHInit(double *var_in1, double *var_in2, double *var_out, int size);

__global__ void SweX(double *Hs_in, double *h, double *Hs_out, double *K2w,
    double *K2e, double *Ztopo, double *mann, double *eRF, double *IN, 
    double *ET, int M, int N, int t);

__global__ void SweY(double *Hs_in, double *h, double *Hs_out, double *K2n, 
    double *K2s, double *Ztopo, double *mann, double *eRF, double*IN, 
    double *ET, int M, int N, int t);

__global__ void SWE_Explicit(double *Hs_in, double *h, double *Hs_out, 
    double *K2w, double *K2e, double *K2n, double *K2s, double *Ztopo, 
    double *mann, double *eRF, double*IN, double *ET, int M, int N, int t);
