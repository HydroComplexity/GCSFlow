#ifndef HOSTDATA2D_H
#define HOSTDATA2D_H
  class HostData2D
  {
    public:
      double *Ztopo;
      double *Hs;
      double *Hs_out;
      double *mann;      
      double *h;
      double *h_out;
      double *IN_h;
      double *Psidiff_h;
      double *u2_h;
      double *v2_h;
      double *Kw_h;
      double *Ke_h;
      double *Kn_h;
      double *Ks_h;
      int *iter_z_h;
      int *iter_x_h;
      int *iter_y_h;
      HostData2D (int sizeZ, int sizeX, int sizeY)
      {
        Ztopo = new double[sizeZ];
        Hs = new double[sizeZ];
        Hs_out = new double[sizeZ];
        mann = new double[sizeZ];
        h = new double[sizeZ];
        h_out = new double[sizeZ];
        IN_h = new double[sizeZ];
        Psidiff_h = new double[sizeZ];
        u2_h = new double[sizeZ];        
        v2_h = new double[sizeZ];
        Kw_h = new double[sizeZ];
        Ke_h = new double[sizeZ];
        Hs = new double[sizeZ];
        Kn_h = new double[sizeZ];
        Ks_h = new double[sizeZ];
        iter_z_h = new int[sizeZ];
        iter_x_h = new int[sizeX];
        iter_y_h = new int[sizeY];
      };
  };
#endif 

#ifndef HOSTDATA3D_H
#define HOSTDATA3D_H
  class HostData3D
  {
    public:
      double *Psi_in;
      double *Psi_out;
      double *theta_out;
      double *Ksat;
      
      HostData3D (int size)
      {
        Psi_in = new double[size];
        Psi_out = new double[size];
        theta_out = new double[size];
        Ksat = new double[size];
      };
  };
#endif 
  

#ifndef HOSTDATATIME_H
#define HOSTDATATIME_H
  class HostDataTime
  {
    public:
      double *PPT_data;
      double *ET_data;
      
      HostDataTime (int size)
      {
        PPT_data = new double[size];
        ET_data = new double[size];
      };
  };
#endif 
  
  
#ifndef HOSTVARTIME_H
#define HOSTVARTIME_H
  class HostVarTime
  {
    public:
      double *PPT;
      double *ET;
      double *IF;
      double *qss;
      
      HostVarTime (int size)
      {
        PPT = new double[size];
        ET = new double[size];
        IF = new double[size];
        qss = new double[size];
      };
  };
#endif 