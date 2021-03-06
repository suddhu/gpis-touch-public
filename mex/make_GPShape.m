EIGEN_PATH = '/usr/include/eigen3/Eigen/src';
SRC = '../cpp/src';
INC = '../cpp/include';

disp('mex mexGPShape.cpp GPShape.cpp...');
mex('mexGPShape.cpp',...
    strcat(SRC,'/GPShape.cpp'),...
    strcat(SRC,'/contours.cpp'),...
    strcat('-I',INC),...
    strcat('-I', EIGEN_PATH),...
    'CXXFLAGS=-std=c++14 -pthread -fPIC',...
    'LDFLAGS="\$LDFLAGS"','-O');
