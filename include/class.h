#include <string>
#include <fstream>

#ifndef PROJECT_H
#define PROJECT_H
  // Project class includes general parameters/options for Project
  class ProjectClass 
  {
    public:
      const char *Name;       // Project name
      double dt_data;         // Time step in data file
      int substeps;           // Number of sub time steps
      int time_data;          
      int time_steps;
      int PrintSubStep;       // Number of time step for printing info
      int SaveStep;           // Number of time step for saving output
      int Save2DOverland;     // Option to save 2D model or not (0 or 1)
      int Save3DSubsurface;   // Option to save 3D model or not (0 or 1)
      int SaveIteration;      // Option to save iteration count or not (0 or 1)
      int SaveInfiltration;   // Option to save infiltration or not (0 or 1)
      int BC2D;
      int BLOCK_DIMX;         // GPU Block dimension in X
      int BLOCK_DIMY;         // GPU Block dimension in X
      int BLOCK_DIMZ;         // GPU Block dimension in X
      int BSZ;                // GPU Block Size
      int TSZ;                // GPU Thread Size
  };
#endif  

#ifndef OVERLANDCLASS_H
#define OVERLANDCLASS_H
  // Overland class includes variables and parameters for overland flow model
  class OverlandClass
  {
    public:
      double Numeric;         // Numeric scheme selection (0: Explicit, 1: ADI)
      double Delta;           // Delta parameters
      double MinimumDepth;    // Minimum depth for wetting and drying
      double CriticalDepth;   // Critical Depth
      double K0;              // K0 parameters
  };
#endif  

#ifndef SUBSURFACECLASS_H
#define SUBSURFACECLASS_H
  // Subsurface class includes variables/parameters for subsurface flow model
  class SubsurfaceClass
  {
    public:
      double Psimin;          // Minimum pressure head
      double stoptol;         // Stop tolerance for convergence
      double PsiBottom;       // Bottom pressure head for boundary
      double PsiNorth;        // North interface pressure head for boundary
      double PsiSouth;        // South interface pressure head for boundary
      double PsiEast;         // East interface pressure head for boundary
      double PsiWest;         // West interface pressure head for boundary
      int BoundBottom;        // Bottom boundary type
      int BoundNorth;         // North interface boundary type
      int BoundSouth;         // South interface boundary type
      int BoundEast;          // East interface boundary type
      int BoundWest;          // West interface boundary type
      int MaxIterX;           // Maximum iteration in X 
      int MaxIterY;           // Maximum iteration in Y 
      int MaxIterZ;           // Maximum iteration in Z 
  };
#endif  

#ifndef SOILPROPERTIESCLASS_H
#define SOILPROPERTIESCLASS_H
  // Soil properties class includes parameters for soil physical properties
  class SoilPropertiesClass
  {
    public:
      double Alpha;           // Alpha in van Genuchten model
      double Theta_R;         // Residual water content
      double Theta_S;         // Saturated water content
      double n;               // Pore-size distribution
      double porosity;        // Soil porosity
      double Ss;              // Specific storage
      double Ksat;            // Saturated hydraulic conductivity
  };
#endif  

#ifndef DOMAINCLASS_H
#define DOMAINCLASS_H
  // Domain class includes info for domain of simulation
  class DomainClass
  {
    public:
      double GridSizeDx;      // Space horizontal resolution
      double GridSizeDy;      // Space horizontal resolution
      double GridSizeDz;      // Space vertical resolution
      int Nx;                 // Number of cell in x-direction
      int My;                 // Number of cell in y-direction
      int Pz;                 // Number of cell in z-direction
  };
#endif  

#ifndef DIMENSION_H
#define DIMENSION_H
  class DimensionClass
  {
    public:
      int dim_topo[2];        // dimensions of topography [M,N]
      int dim_para[2];        // dimensions of parameters [M,N]
      int dim_Psi[3];         // dimensions of pressure head variable [M,N,P]
      int dim_Ksat[3];        // dimensions of pressure head variable [M,N,P]
      int dim_forcing[1];     // dimensions of forcing [Time]
  };
#endif  
  
#ifndef FILENAMECLASS_H
#define FILENAMECLASS_H
  class FileNameClass
  {
    public:
      const char *topo_file;          // Topography filename
      const char *init_conds_file;    // Initial condition filename  
      const char *forcing_file;       // Forcings filename
      const char *parameter_file;     // Parameters filename  
      const char *output3D;           // Output filename for 3D variables
      const char *output2D;           // Output filename for 2D variables
      const char *output_veloc2D;     // Output filename for 2D velocity
      const char *output_K2D;
      const char *output_h;
      const char *output_u;
      const char *output_v;
      const char *output_qss;
  };
#endif
  