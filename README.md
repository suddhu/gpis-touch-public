# GP implicit surfaces for shape estimation

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) &nbsp; <img height="20" src="https://rpl.ri.cmu.edu/images/rpl4_cropped.png" alt="RPL-logo" />

![butter_2d](/media/butter_2D.gif) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ![butter_3d](/media/butter_3D.gif)

This library generates an implicit surface representation from sparse tactile information, used for the paper [Tactile SLAM: Real-time inference of shape and pose from planar pushing
](https://arxiv.org/pdf/2011.07044.pdf)

## Folders
- **cpp/**: The C++ .cpp and .h files
- **data/**: Has contact and shape data
  - `python get_contact.py --shape rect1` lets you select contact points and normals to make your own dataset. First click sets the contact point, second click draws the normal. Normal must be *away* from the object.
  - **shapes/**: the ground truth shape `.mat` files and scripts
  - **contacts/**: The generated contact files are stored here
- **matlab/**:  `GPshape.m` is the launch file and `viz_shape.m` has the visualization code
  - **standalone/**: old MATLAB code which implements a basic version of GPIS
- **mex/**: Contains the makefile, mexfile, and `test_gp.cpp`

## CMakeLists Notes 
- Modify the `EIGEN3_INCLUDE_DIR` in `mex/CMakeLists.txt` and `cpp/CMakeLists.txt`
- Modify `EIGEN_PATH` in `mex/make_GPShape.m`
- Modify `BOOST_ROOT` in `cpp/CMakeLists.txt`
## Mex executable
Compile and run from MATLAB
### Compile
```
cd mex/
make_GPShape
```
### Run
```
cd matlab/
GPshape
```

## C++ executable 
### Compile
```
cd mex/
mkdir build
cd build/
cmake ..
make -j
```
### Run
```
./test_gp
```

## Reference and Acknowledgements
- The mex functions have been adapted from Lee, Bhoram, et al. "Online continuous mapping using gaussian process implicit surfaces." 2019 International Conference on Robotics and Automation (ICRA). IEEE, 2019 ([github](https://github.com/leebhoram/GPisMap))
- Another useful open-source reference was "Dexterous grasping under shape uncertainty", Miao Li, Kaiyu Hang, Danica Kragic and Aude Billard, Robots and Autonomous Systems, 2015 ([github](https://github.com/epfl-lasa/GPIS))
- The contouring functions are from [conrec](http://paulbourke.net/papers/conrec/).
- C++ library to generate mesh grid: [meshgen](https://github.com/xiaohongchen1991/meshgen)
- Shape models taken from the [MIT push dataset](https://mcube.mit.edu/push-dataset/index.html)
- cpp plotting uses [matplotlib-cpp](https://github.com/lava/matplotlib-cpp) 
- `kdtree_eigen.h` is developed by Fabian Meyer based on "Analysis of Approximate Nearest Neighbor Searching with Clustered Point Sets" by Songrit Maneewongvatana and  David M. Mount

## Citation 
Feel free to use the library as you please. If you find it helpful, please consider referencing: 

```BibTeX
@article{suresh2020tactile,
  title={Tactile SLAM: Real-time inference of shape and pose from planar pushing},
  author={Suresh, Sudharshan and Bauza, Maria and Yu, Kuan-Ting and Mangelson, Joshua G and Rodriguez, Alberto and Kaess, Michael},
  journal={arXiv preprint arXiv:2011.07044},
  year={2020}
}
```

## Roadmap
- [ ] 3D object reconstruction 
- [ ] More kernels 
