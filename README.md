GCSFlow Model v1.0
=======
GPU-based Conjunctive Surface Sub-surface Flow Model

[![Build Status](https://travis-ci.org/simkimsia/UtilityBehaviors.png)](https://travis-ci.org/simkimsia/UtilityBehaviors)


This is an open source, GPU-based program for modeling conjunctive 2D surface and 3D sub-surface flow. The model is written in CUDA C++ and has been tested under GNU/Linux.

### Documentation

* Under review for publication


### Installation
```
make build
make
```
<small>See or edit Makefile to get or change the executable file for the GCSFlow model  
Note: NetCDF library and GPUs with computability 2.0 or greater are required to compile the code.</small>



### Run  
```
cd Example_folder
./executable ../path_to_config/config_file.cfg
```


### License
This software is freeware and is released under LGPL. see LICENSE file for more information. 


### Contact Authors
* Phong Le: <mailto:levuvietphong@gmail.com>
* Praveen Kumar: <mailto:kumar1@illinois.edu>  
Questions and suggestions are welcome.