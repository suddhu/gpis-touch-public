/*
 * get_shapes.hpp
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
#include <cmath>
#include <vector>
#include <unordered_map>

#ifndef __GET_SHAPES_HPP__
#define __GET_SHAPES_HPP__

using namespace std;

void ground_truth_shape(vector<vector<double>>& shape_data, string shape_id);
// https://gist.github.com/damithsj/c96a8482b282a3dc89bd
vector<double> linspace(double min, double max, int n);
vector<vector<double>> makeShapePolyRect(double longSide, double shortSide);
vector<vector<double>> makeShapePolyTri(double shortSide1, double shortSide2, double longSide);
vector<vector<double>> makeShapeEllip(double a, double b);
vector<vector<double>> makeShapePolyNGon(double side, int n);
#endif
