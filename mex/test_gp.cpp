/*
 * test_gp.cpp: cpp test code for GPShape library
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

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "GPShape.hpp"
#include "get_shapes.hpp"
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

// plot the intermediate GP results with contact data
void plotResult(cnt contour, cnt ground_truth, Eigen::Vector2d contact_point, Eigen::Vector2d contact_normal, string shape_id, int ts) {

    bool pause = false;
    vector<double> shape_x, shape_y;
    vector<double> ground_x, ground_y;
    vector<double> c_x, c_y, n_x, n_y; // for plotting

    for (int i = 0; i < contour.size(); i++) {
        shape_x.push_back(contour[i][0]); 
        shape_y.push_back(contour[i][1]); 
    }
    shape_x.push_back(shape_x[0]); 
    shape_y.push_back(shape_y[0]); 

    for (int i = 0; i < ground_truth.size(); i++) {
        ground_x.push_back(ground_truth[i][0]); 
        ground_y.push_back(ground_truth[i][1]); 
    }

    ground_x.push_back(ground_x[0]);
    ground_y.push_back(ground_y[0]);

    // for plotting
    c_x.push_back(contact_point(0)); c_y.push_back(contact_point(1)); 
    n_x.push_back(contact_normal(0)); n_y.push_back(contact_normal(1)); 

    // TODO: plot wrt object pose 
    plt::plot(shape_x, shape_y, "g-"); // optimized shape

    plt::plot(ground_x, ground_y, "k--"); // ground truth shape
    plt::plot(c_x, c_y, "bx"); // contact points (all)

    plt::quiver(c_x, c_y, n_x, n_y);

    string title_string = shape_id + string(" t = ") + to_string(ts); 
    plt::title({title_string});
    plt::axis("equal"); 

    plt::xlim(-0.09, 0.09); 
    plt::ylim(-0.09, 0.09);

    if(pause){
        plt::show();
    }
    else {     
        plt::show(false);
        plt::pause(0.0001);
        plt::clf();
    }
}

int main(int argc, char** argv){

    // ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-rect1-20200115-1026.txt");
    ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-rect1-20200810-1811.txt");
    // ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-rect1-20200810-1616.txt");
    // ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-rect1-20200831-1244.txt");
    // ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-rect1-20200831-1248.txt");
    // ifstream infile("/home/suddhu/software/GPIS/data/contacts/contacts-ellip2-20200506-2345.txt");

    string shape_id = "rect1";
    // string shape_id = "ellip2";

    cnt ground_truth; 
    ground_truth_shape(ground_truth, shape_id); 

    // initialize GP
    std::vector<double> pvar = {1e-5, 1e-5, 1e-5}; // meas variance: scalar
    std::vector<double> priorNoise = {1e-2, 1e-2, 1e-2}; // meas variance: scalar
    double testLim = 0.1; double testRes = testLim*0.1; double priorRad = 0.03;
    int gp_count = 16; // 4*4
    std::string kernel = "thinplate";    // kernel options: "thinplate", "gaussian", "matern"

    std::shared_ptr<GPShape> gpshape; // shape object
    gpshape = std::make_shared<GPShape>(pvar, priorNoise, kernel, testLim, testRes, gp_count, true);
    // GPShape(std::vector<double> var_, std::vector<double> prior_var_, std::string kernel_, const double& test_lim_, const double& res_, const int& gp_count_) : 

    // add circular prior to all GPs
    gpshape->addPrior(priorRad, 4, Eigen::Vector3d::Zero());

    string line_data;
    cnt contour; 
    Emx fMean, fVar; Evx contourVar;

    int t = 0;
    while (getline(infile, line_data))
    {
        stringstream line_stream(line_data);
        double px, py, nx, ny; char comma;
        if (line_stream >> px >> comma >> py >> comma >> nx >> comma >> ny) {
            Eigen::Vector2d contact_point(px, py);
            Eigen::Vector2d contact_normal(nx, ny);

            bool updated = gpshape->update(contact_point, contact_normal, Eigen::Vector3d::Zero());

            if (updated)    gpshape->test(contour, fMean, fVar, contourVar);
            else            std::cout<<"skip testing!\n";
            plotResult(contour, ground_truth, contact_point, contact_normal, shape_id, t);
        }
        t++;
    }
    // gpshape->printTiming();
    return 0;
}
