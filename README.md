<img src="https://cloud.githubusercontent.com/assets/1987216/4551081/ee29a83a-4e6c-11e4-8707-355641753c64.png" width="500px" hspace=5>  
GPU-based Conjunctive Surface Sub-surface Flow Model v1.0

[![Build Status](https://travis-ci.org/simkimsia/UtilityBehaviors.png)](https://travis-ci.org/simkimsia/UtilityBehaviors)


This is an open source, GPU-based program for modeling conjunctive 2D surface and 3D sub-surface flow. The model is written in CUDA C++ and has been tested under GNU/Linux.

### Documentation
Under review for publication


### Installation
Prerequisites
```
NetCDF 4.0 or newer
CUDA 5.0 or newer
g++ 4.2 or newer
```

To build
```
$ make build
$ make
```


### Run  
```
$ cd bin
$ ./gcsflow ../Test_Cases/path_to_test/config_filename.cfg
```
Edit Makefile to change the compiling options.  
Model options and parameters are in configuration files.

### License
This software is freeware and is released under LGPL. See LICENSE file for more information. 


### Contact Authors
* Phong Le: <mailto:levuvietphong@gmail.com>
* Praveen Kumar: <mailto:kumar1@illinois.edu>

Questions and suggestions are welcome.