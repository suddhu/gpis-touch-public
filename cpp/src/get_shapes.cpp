/*
 * get_shapes.cpp: Get shape models for test_gp.cpp
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

#include "get_shapes.hpp"
using namespace std;

std::unordered_map<std::string, std::initializer_list<double>> shapeid2dim = {
        {"rect1",{0.090, 0.090}},
        {"rect2",{0.08991, 0.11258}},
        {"rect3",{0.13501, 0.08994}},
        {"tri1",{0.12587, 0.12590, 0.178}},
        {"tri2",{0.12587, 0.15100, 0.1962}},
        {"tri3",{0.12561, 0.1765, 0.2152}},
        {"ellip1",{0.105, 0.105}},
        {"ellip2",{0.105, 0.13089}},
        {"ellip3",{0.105, 0.157}},
        {"hex",{0.06050, 6.0}}
    };
    
void ground_truth_shape(vector<vector<double>>& shape_data, string shape_id){

    std::vector<double> dims;

    if((shape_id == "rect1") || (shape_id == "rect2") || (shape_id == "rect3")) {
      for (auto d : shapeid2dim[shape_id])
        dims.push_back(d);
      shape_data = makeShapePolyRect(dims[0], dims[1]);
    } 
    else if ((shape_id == "tri1") || (shape_id == "tri2") || (shape_id == "tri3")) {
      for (auto d : shapeid2dim[shape_id])
        dims.push_back(d);
      
      shape_data = makeShapePolyTri(dims[0], dims[1], dims[2]);
    }
    else if ((shape_id == "ellip1") || (shape_id == "ellip2") || (shape_id == "ellip3")) {
      for (auto d : shapeid2dim[shape_id])
        dims.push_back(d);

      shape_data = makeShapeEllip(dims[0], dims[1]);
    }
    else if (shape_id == "hex") {
      for (auto d : shapeid2dim[shape_id])
        dims.push_back(d);
      
      shape_data = makeShapePolyNGon(dims[0], (int)dims[1]);
    }
}

// https://gist.github.com/damithsj/c96a8482b282a3dc89bd
vector<double> linspace(double min, double max, int n)
{
	vector<double> result;
	// vector iterator
	int iterator = 0;

	for (int i = 0; i <= n-2; i++)	
	{
		double temp = min + i*(max-min)/(floor((double)n) - 1);
		result.insert(result.begin() + iterator, temp);
		iterator += 1;
	}

	//iterator += 1;

	result.insert(result.begin() + iterator, max);
	return result;
}

vector<vector<double>> makeShapePolyRect(double longSide, double shortSide) {
    double a = longSide / 2.0;
    double b = shortSide / 2.0;

    vector<vector<double>> polyrect; 
    polyrect.push_back(vector<double>{a,b}); 
    polyrect.push_back(vector<double>{-a,b}); 
    polyrect.push_back(vector<double>{-a,-b}); 
    polyrect.push_back(vector<double>{a,-b}); 
    return polyrect;
}

vector<vector<double>> makeShapePolyTri(double shortSide1, double shortSide2, double longSide){
    double a = shortSide1;
    double b = shortSide2;
    double c = longSide;
    double d = 0.090 / 2.0;  // from the rectangle coordinate system
    
    vector<vector<double>> polytri; 
    polytri.push_back(vector<double>{d,d}); 
    polytri.push_back(vector<double>{d-b,d}); 
    polytri.push_back(vector<double>{d,d-a}); 
    return polytri;
}

vector<vector<double>> makeShapeEllip(double a, double b){
    vector<double> gamma = linspace(0, 2*M_PI, 100);

    double ax1 = a/2.0; 
    double ax2 = b/2.0; 
    
    vector<vector<double>> ellip; 
    for (int i = 0; i < gamma.size(); i++) {
        ellip.push_back(vector<double>{ax1*cos(gamma[i]), ax2*sin(gamma[i])});
    }

    return ellip;
}

vector<vector<double>> makeShapePolyNGon(double side, int n){
    vector<vector<double>> polyngon; 
    for (int i = 0; i < n; i++) {
        double theta = (2*M_PI/n)*i;
        polyngon.push_back(vector<double>{side*cos(theta), side*sin(theta)});
    }
    return polyngon;
}
