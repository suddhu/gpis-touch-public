/*
 * mexGPShape.cpp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License v3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of any FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License v3 for more details.
 *
 * You should have received a copy of the GNU General Public License v3
 * along with this program; if not, you can access it online at
 * http://www.gnu.org/licenses/gpl-3.0.html.
 *
 * Author: Sudharshan Suresh <suddhu@cmu.edu>
 */

/* 
The mex interface code was adapted and repurposed from the GPisMap codebase (https://github.com/leebhoram/GPisMap): 
Lee, Bhoram, et al. "Online continuous mapping using gaussian process implicit surfaces." 2019 International Conference on Robotics and Automation (ICRA). IEEE, 2019.
 */

#include "mex.h"
#include <vector>
#include <memory>
#include <iostream>
#include <chrono>
#include "GPShape.hpp"

static GPShape* gp = 0;

static double t_dur = 0.0;
static std::chrono::high_resolution_clock::time_point t_start;
static std::chrono::high_resolution_clock::time_point t_end;

template <class T, class I = unsigned>
struct MatrixWrapper
{
    T *data;
    I nrows, ncols;

    MatrixWrapper() 
        { clear(); }
    MatrixWrapper( T *data_, I nrows_, I ncols_ )
        { set(data_,nrows_,ncols_); }

    inline void clear()
        { data = NULL; nrows = ncols = 0; }

    inline void set( T *data_, I nrows_, I ncols_ )
        { data = data_; nrows = nrows_; ncols = ncols_; }

    inline T& operator() ( I r, I c ) const
        { return data[ r + nrows*c ]; }
};

mxArray * getMexArray(const std::vector<double>& v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

// https://www.mathworks.com/help/matlab/apiref/mexfunction.html?searchHighlight=mexfunction
// nlhs: # of output arguments
// plhs: Array of pointers to the expected mxArray output arguments.
// nrhs: # of input arguments
// prhs: Array of pointers to the mxArray input arguments. 
void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {

    // get flag
    char command[128];
    mxGetString(prhs[0],command,128);
    std::string commstr(command);

    if (commstr.compare("init")==0) {

        if (gp !=0) {
            delete gp;
            gp = 0;
            std::cout<<"init(): GP already intialized! Clearing..."<< std::endl;
        }
        if ((nrhs == 7) && (gp == 0) ){ 
            double *pvar, *priorvar;
            pvar = mxGetPr(prhs[1]); // meas variance: scalar
            priorvar = mxGetPr(prhs[2]); // prior variance: scalar
            double testLim = (double)mxGetScalar(prhs[3]); // kernel param: scalar
            double testRes = (double)mxGetScalar(prhs[4]); // kernel param: scalar
            double priorRad = (double)mxGetScalar(prhs[6]); // kernel param: scalar

            bool isLocal = false; int GPcount = 1;
            if ((int)mxGetScalar(prhs[5]) !=0) {
                isLocal = true;
                GPcount = 16;
            }

            std::string kernel = "thinplate";
            gp = new GPShape(std::vector<double>(pvar, pvar + 3), std::vector<double>(priorvar, priorvar + 3), kernel, testLim, testRes, GPcount, isLocal);
            gp->addPrior(priorRad, 4, Eigen::Vector3d::Zero());
            return;
        } 
        else {
            std::cout<<"init(): incorrect parameters"<< std::endl;
            return;     
        }
    }

    else if (commstr.compare("update")==0) {
        if (nrhs == 3)
        {
            double * px = (double *)mxGetData(prhs[1]); // training input (x,y): 2 x 1
            double * py = (double *)mxGetData(prhs[2]); // training output (nx, ny): 2 x 1
            t_start = std::chrono::high_resolution_clock::now();

            Eigen::Vector2d contact_point(px[0], px[1]);
            Eigen::Vector2d contact_normal(py[0], py[1]);

            gp->update(contact_point, contact_normal, Eigen::Vector3d::Zero());
            
            t_end= std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t_end-t_start); // reset
            t_dur = time_diff.count();
            // plhs[0] = mxCreateLogicalScalar(updated);	
        } 
        else if (gp ==0) {
            std::cout<<"update(): GP not intialized!"<< std::endl;
            return;
        } 
        else {
            std::cout<<"update(): incorrect parameters"<< std::endl;
            return;     
        }
    }
    else if (commstr.compare("test")==0) {
        if ((gp != 0) && (nrhs == 1)){
            t_start = std::chrono::high_resolution_clock::now();

            cnt contour; 
            Emx fMean, fVar; Evx contourVar;
            gp->test(contour, fMean, fVar, contourVar);
            t_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t_end-t_start); // reset
            t_dur = time_diff.count();
            
            // return contour
            int nrows = contour.size(); 
            int ncols = 2; 
            double *M_data = new double[ nrows*ncols ];             // allocate a temporary matrix
            MatrixWrapper<double> M( M_data, nrows, ncols );

            for(int j = 0; j<nrows; j++) {
                M(j, 0) = (double) contour[j][0];
                M(j, 1) = (double) contour[j][1];
            }             // access values with M(i,j) as a "normal" matrix

            plhs[0]  =  mxCreateDoubleMatrix( nrows, ncols, mxREAL );
            MatrixWrapper<double> out0( mxGetPr(plhs[0]), nrows, ncols );
            for ( unsigned c = 0; c < ncols; c++ ) {
                for ( unsigned r = 0; r < nrows; r++ )
                    out0(r,c) = M(r,c);
            }             // copy to the Mex output

            M.clear();             // free temporary allocation
            delete[] M_data;

            // return fMean transposed
            mwSize mean_rows = fMean.rows();
            mwSize mean_cols = fMean.cols();
            plhs[1] = mxCreateDoubleMatrix(mean_rows, mean_cols, mxREAL); // Create MATLAB array of same size
            Eigen::Map<Eigen::MatrixXd> map1(mxGetPr(plhs[1]), mean_rows, mean_cols); // Map the array
            map1 = fMean; // Copy

            // return fVar transposed
            mwSize var_rows = fVar.rows();
            mwSize var_cols = fVar.cols();
            plhs[2] = mxCreateDoubleMatrix(var_rows, var_cols, mxREAL); // Create MATLAB array of same size
            Eigen::Map<Eigen::MatrixXd> map2(mxGetPr(plhs[2]), var_rows, var_cols); // Map the array
            map2 = fVar; // Copy
        }
        else if (gp==0) {
            std::cout<<"test(): GP not intialized!"<< std::endl;
            return;
        }
        else {
            std::cout<<"test(): incorrect parameters"<< std::endl;
            return;     
        }
    }
    else if (commstr.compare("reset")==0) {
        if (gp != 0){
            delete gp;
            gp = 0;
            std::cout << "GP cleared\n" << std::endl;
        }
    }
    return;
}
