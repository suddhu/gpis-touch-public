/*
 * GPShape.hpp
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
 * Author: Sudharshan Suresh <suddhu@cmu.edu>, 2020
 */

#ifndef __GP_SHAPE_HPP__
#define __GP_SHAPE_HPP__

#include <vector>
#include <deque>
#include <memory>
#include <iostream>
#include <cstdint>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <math.h>
#include <functional>
#include <numeric>
#include <chrono>
#include <thread>
#include <algorithm>

#include "meshgrid.hpp"
#include "contours.hpp"
#include "kdtree_eigen.h"

typedef Eigen::MatrixXd Emx;
typedef Eigen::VectorXd Evx;


// https://www.learncpp.com/cpp-tutorial/8-16-timing-your-code/
class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;
	std::chrono::time_point<clock_t> m_beg;
 
public:
	Timer() : m_beg(clock_t::now()){}
	void reset() {m_beg = clock_t::now();}
	double elapsed() const {
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

// structure stores x, y contact positions and nx, ny normal directions
struct Contact
{
    double x, y, nx, ny;

    Contact(double _x, double _y, double _nx, double _ny): x(_x), 
    y(_y), nx(_nx), ny(_ny){}
    Contact():  x(0.0), y(0.0), nx(0.0), ny(0.0){}

    Contact & operator = (const Contact & that) {
        x = that.x; y = that.y;
        nx = that.nx; ny = that.ny;
        return *this;
    }

    std::ostream& operator<<(std::ostream& stream) {
        std::cout<<"x: "<<x<<" y: "<<y<<" nx: "<<nx<<" ny: "<<ny<<std::endl;
        return stream;
    }
};

// individual GP structure
struct GP
{
    std::vector<Contact> measurements; // measurements 
    Emx outputs = Emx::Zero(0, 1); //outputs

    Evx currPose; // pose to transform 
    int numMeasurements, numPrior; // counters
    Emx Kxx, Kxx_L, alpha;         // Kernels, cholesky 

    std::vector<double> gp_center;
    int gp_count, gp_id;

    GP() {}
    GP(int gp_id_, int gp_count_, double test_lim_): numMeasurements(0), numPrior(0),
    gp_count(gp_count_), gp_id(gp_id_) {computeGPCenter(test_lim_);}
    ~GP(){};

    void computeGPCenter(double test_lim_) {
        int length = (int)sqrt(gp_count); 
        double dist = 2*test_lim_/((double)length + 1.0);
        int r = gp_id/length; int c = gp_id % length;
        gp_center = {-test_lim_ + (double)(r + 1)*dist, -test_lim_ + (double)(c + 1)*dist};
    };

    std::vector<double> getGPCenter() { return gp_center;};
};

class GPShape{

private:
    std::vector<double> var, prior_var; // variance
    double vThresh; //accept thresh

    std::string kernel;     // choice of kernel function
    double R;               // Thin plate kernel parameter
    double l = 0.25;        // Gaussian kernel parameter
    double rho = 0.15;      // Matern kernel parameter
    double sigma = 0.01;    // Kernel parameter

    double gp_radius; // local GP radius
    int gp_count; // total # of GPs

    double testLim, res; // grid limit, res
    meshgen::mesh_grid<double, 0, 2> x_test;
    meshgen::mesh_grid<double, 1, 2> y_test;
    double* xx; double* yy;
    int test_side2, test_side; // eg: 4 x 4 grid
    
    Emx M_test, C_test; 
    std::vector<GP*> gps; // GP class object eg: gps[5] = (2, 1)
    std::vector<std::vector<double>> gps_centers; 

    kdt::KDTreed::MatrixI nnGP; // nearest neighbor GP
    std::vector <std::vector <int> > nnGPValid; //valid nearest neighbor

    bool noPreFilter = false;
    bool localExperts;
    std::vector<bool> isUpdated;
    std::vector<double> update_time, test_time;

public:

    GPShape(std::vector<double> var_, std::vector<double> prior_var_, std::string kernel_, const double& test_lim_, const double& res_, const int& gp_count_, bool localExperts_ = true) : 
    var(var_), prior_var(prior_var_),
    vThresh(1e-1*var_[0]), kernel(kernel_),
    gp_count(gp_count_), testLim(test_lim_), res(res_), localExperts(localExperts_) {    

        // R = 2*testLim;
        if (localExperts)   R = 1.5*testLim;
        else                R = 10*testLim; // R must be greater than the greatest distance between measurements

        double x_start = -test_lim_; double y_start = -test_lim_;
        double x_end = test_lim_; double y_end = test_lim_;

        int N_x = (int)((x_end - x_start)/res);
        int N_y = (int)((y_end - y_start)/res);
        std::vector<double> x(N_x + 1); meshgen::linspace(x.begin(), x.end(), x_start, res);
        std::vector<double> y(N_y + 1); meshgen::linspace(y.begin(), y.end(), y_start, res);
        std::tie(x_test, y_test) = meshgen::meshgrid(x.begin(), x.end(), y.begin(), y.end()); 

        test_side = static_cast<int>(x_test.size1()); // number of test elements 
        test_side2 = test_side*test_side; 

        M_test  = Emx::Zero (test_side, test_side); 
        C_test  = Emx::Zero (test_side, test_side); 

        xx = new double[test_side];
        yy = new double[test_side];

        for (int i = 0; i < test_side; ++i) { xx[i] = x_test(i, 0); yy[i] = y_test(0, i);}

        initGPs();  buildAndQueryTree();  

        Eigen::initParallel();
    }

    ~GPShape(){};

    void buildAndQueryTree();
    void initGPs();

    // update kernels and params 
    void updateKxx(const int& idx);
    void updateAlpha(const int& idx);

    void getKxStar2(const Emx& x_star, Emx& K_xStar_xStar) const;
    void getKxStar2(const vector<double>& x, Emx& K_xStar_xStar) const;

    void getKxStar(const int& idx, const Emx& x_star,  Emx& KxStar) const;
    void getKxStar(const int& idx, const vector<double>& x,  Emx& KxStar) const;


    void addKernelNoise(Emx& kernel, const int& num_prior);
    void getKxCheck(const int& idx, const Emx& x_check, Emx& KuStar) const;

    std::vector<std::vector<double>> getCenters() const {return gps_centers;}; 

    bool proximityCheck(const int& idx, const Evx& x) const;

    Emx kernelBlock(const std::vector<double>& x1, const std::vector<double>& x2) const;
    Emx kernelBlockThinPlate(const std::vector<double>& x1, const std::vector<double>& x2) const;
    Emx kernelBlockGaussian(const std::vector<double>& x1, const std::vector<double>& x2) const;
    Emx kernelBlockMatern(const std::vector<double>& x1, const std::vector<double>& x2) const;

    // cholesky solutions
    Emx solveCholesky(const Emx& A, const Emx& b) const;
    Emx solveCholeskyWithL(const Emx& L, const Emx& b) const;
    Emx getCholeskyL(const Emx& A) const;
    void updateCholeskyL(const int& idx, const Emx& A);

    // data 
    void addMeasurement(const int& idx, const Evx& c,  const Evx& n);

    void printSize(const Emx& M, const std::string& M_name) const {std::cout<<M_name<<".size() : "<<M.rows()<<" * "<<M.cols()<<std::endl;};

    void printTiming() const {
        std::cout<<"\n\nupdate time:"<<std::endl;
        for(uint i = 0; i < update_time.size(); i++)
            std::cout<<update_time[i]<<std::endl;
        std::cout<<"\n\ntest time:"<<std::endl;
        for(uint i = 0; i < test_time.size(); i++)
            std::cout<<test_time[i]<<std::endl;
    };

    // GP 
    bool preFilter(const int& idx, const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p); 
    void updateThread(const int& start_id, const int& end_id,  const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p, std::vector<bool>& isUpdate);
    bool update(const Evx& contact_pt,  const Evx& contact_normal, const Evx& p);
    bool updateLocal(const int& idx, const Evx& contact_pt,  const Evx& contact_normal, const Evx& p, const bool& isPrior = false);
    void addPrior(const double& r, const int& n, const Evx& p);

    void test(cnt& contour, Emx& m_test, Emx& c_test, Evx& contourVar);
    void testThread(const int& thread_id, const int& start_id,
                     const int& end_id, const int& sz, Emx& m_, Emx& c_) const;
    void testMeasurement(const int& idx, const Evx& p, const Evx& n, 
                    double& v, double& n_dot, double& d) const;
    Evx getContourVariance(const cnt& contour, const Emx& c_test) const;
 
    // misc
    void transformContour(cnts& contours,  const Evx& p) const;
    Evx transformPointFrom(const Evx& pt, const Evx& pose2) const;    // From local to global frame
    Evx transformPointTo(const Evx& pt, const Evx& pose2) const;    // From global to local frame
    double Norm( const std::vector<double>& x,  const std::vector<double>& y) const;    // Euclidean distance norm Norm(x - y)
};


#endif