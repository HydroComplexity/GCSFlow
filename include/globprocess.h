__global__ void CopyVariable(double *var_in, double *var_out, int size);

__global__ void TopForcing(double ppt, double *eff_rain, int size);

__global__ void PondHeadInit(double *ph, int size);

__global__ void TopBottomBound2D(double *Hs, double *Ztopo, double *K2n, 
    double *K2s, int BC2D, int M, int N);

__global__ void LeftRightBound2D(double *Hs, double *Ztopo, double *K2e,
    double *K2w, int BC2D, int M, int N);
 
__global__ void EstimateVelocity(double *u, double *v, double *K2w, double *K2e,
    double *K2n, double *K2s, double *Hs, double *h, int M, int N);

__global__ void GetOutlet(double *h, double *houtlet, double *u, double *uout, 
    double *v, double *vout, int M, int N, int t);

__global__ void getqss(double *IN, double *qss, int N, int t);

__global__ void VarPrint(double *Var, int M, int N, int P);