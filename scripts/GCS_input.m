%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%           GPU-based Conjunctive Surface-Subsurfae Flow Model          %%
%%             Processing topography and forcing data                    %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
% This matlab program is used to provide the following files              %
% in NetCDF format for the GCS-Flow Model:                                % 
%   - Topography: topography.nc                                           % 
%       + Includes the ground elevation surface                           %
%                                                                         %
%   - Initial condition: init_conds.nc                                    %
%       + Includes the 3D initial pressure head condition in the domain   %  
%         of computation.                                                 %
%                                                                         %
%   - Parameters: parameters.nc                                    %
%       + Includes the 3D initial pressure head condition in the domain   %  
%         of computation.                                                 %
%                                                                         %
%   - Forcing data: forcing.nc                                            %
%       + Includes Precpitation and Evaporation over time                 %  
%-------------------------------------------------------------------------%
%   Created by  : Phong Le                                                %
%   Date        : June 09, 2014                                           %
%-------------------------------------------------------------------------%
%=====================================================
% LOCAL VARIABLES.................
%   M: Number of grid cells in X-direction [-]
%   N: Number of grid cells in Y-direction [-]
%   P: Number of grid cells in Z-direction [-] (depth)
%   Ztopo: 2D Ground elevation [L]
%   Psi_init: 3D Initial pressure head [L]
%   PPT_in: 1D Precipitation data [L/T]
%   Evap_in: 1D Evaporation data [L/T]
%=====================================================

  clear all       % Clear all memory and variables      
  close all     	% Close all programs and functions    
  clc

%---------------------------------------------------------
% CHANGE THE CODE HERE TO CREATE APPROPRIATE INPUT FILES
%---------------------------------------------------------
  
  %. . . There are several ways to load data into matlab. 
  % For TOPOGRAPHY, we need to create/load a variable named Ztopo in 
  % matlab workspace. This scripts will write this variable into a 
  % NetCDF file (topography.nc) for GCS-Flow.
  
  % For FORCING, we create/load two variables named PPT_in and Evap_in
  % in matlab workspace. This scripts will write these variables into a 
  % NetCDF file (forcing.nc) for GCS-Flow.

  % For INITIAL CONDITIONS, we create/load a variable named Psi_init
  % in matlab workspace. This scripts will write this variable into a 
  % NetCDF file (init_conds.nc) for GCS-Flow.
  
  % We provide three examples below for reference. Uncomment the method 
  % that works for you to use. You may need to write your own code to 
  % load your data in other formats
  
  %----------------------------------
  %...Example 1: Load from *.mat file
    % This *.mat_file may include also data needed. If not, continue loading
    % other files for other variables.
  
    %load mat_file
    %load topo_30
    %Ztopo = ncread('simple_xy.nc','data');
    
  %----------------------------------
  %...Example 2: Load from *.tif file
    % This tif_file usually includes topograhy data (with georeference).
    % We may need to load other formats for forcing data.
  
    [Z, R, bbox] = geotiffread('tiff/box3.tif');
    Ztopo = double(Z(2:41,2:51));

  %----------------------------------  
  %...Example 3: Load from another NetCDF file
    % We may load multi variables in one NetCDF files.
        
    %Ztopo = ncread('filename.nc','Topo_variable');
    %PPT_in = ncread('filename.nc','PPT_variables');
  
  %----------------------------------  
  % Number of vertical cells. 
  P = 10;
  
  % Get the size of the domain in X & Y directions  
  [M, N] = size(Ztopo);              
    
  %%...Pressure head initialization . . .
  Psi_init = zeros(M,N,P);            % Initialize the pressure head
  Psi_vert = linspace(-1.0, -1.0, P);
  for j=1:M
    for i=1:N
      Psi_init(j,i,:) = Psi_vert;
    end
  end

  %%...Parameters . . .
  mann = 0.25*ones(M,N);

  %%...Forcing data . . .  
  num_timesteps = 1000;
  load PPT.dat
  PPT_in = PPT(1:num_timesteps)/1000;
  PPT_sum 	= sum(PPT_in) * 1000
  %PPT = 0.33*ones(num_timesteps,1);
  %PPT(201:end) = 0;
  

  Evap_in = 0.0000*ones(num_timesteps,1);
  
  DIR = '../data/';
  %...Check if RESULTS folder exits
  if (exist(DIR,'dir') ~= 7)
    mkdir(DIR);
  end
  
  topo_file = 'topography.nc';
  init_conds_file = 'init_conds.nc';
  forcing_file = 'forcings.nc';
  parameter_file = 'parameters.nc';
  
%---------------------------------------------------------
% DO NOT CHANGE THE CODE BELOW THIS LINE 
%---------------------------------------------------------
 
  %%... Creating topography data
  if (exist([DIR,topo_file], 'file') == 2)    % If file exist, remove file
    delete([DIR,topo_file])
  end
  nccreate([DIR,topo_file], 'Ztopo','Dimensions', {'x' N 'y' M}, 'Datatype' , 'double');
  ncwrite([DIR,topo_file],'Ztopo', Ztopo');
  
  %%... Creating initial condition file
  if (exist([DIR,init_conds_file], 'file') == 2)    % If file exist, remove file
    delete([DIR,init_conds_file])
  end
  nccreate([DIR,init_conds_file], 'Psi_init','Dimensions', {'x' N 'y' M 'z' P}, 'Datatype' , 'double');
  ncwrite([DIR,init_conds_file],'Psi_init', permute(Psi_init,[2 1 3]));
  
    %%... Creating topography data
  if (exist([DIR,forcing_file], 'file') == 2)    % If file exist, remove file
    delete([DIR,forcing_file])
  end
  nccreate([DIR,forcing_file], 'PPT_in','Dimensions', {'row' num_timesteps}, 'Datatype' , 'double');
  nccreate([DIR,forcing_file], 'Evap_in','Dimensions', {'row' num_timesteps}, 'Datatype' , 'double');
  nccreate([DIR,forcing_file], 'num_steps','Datatype' , 'int32');
  
  ncwrite([DIR,forcing_file],'PPT_in', PPT_in);
  ncwrite([DIR,forcing_file],'Evap_in', Evap_in);
  ncwrite([DIR,forcing_file],'num_steps', num_timesteps);  
  
  %%... Creating parameter files
  if (exist([DIR,parameter_file], 'file') == 2)    % If file exist, remove file
    delete([DIR,parameter_file])
  end
  nccreate([DIR,parameter_file], 'mann','Dimensions', {'x' N 'y' M}, 'Datatype' , 'double');
  ncwrite([DIR,parameter_file],'mann', mann');