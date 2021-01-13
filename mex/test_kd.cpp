/*
 * test_kd.cpp: Test code for KDTree library
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
#include "kdtree_eigen.h"

typedef Eigen::MatrixXd Matrix;

int main()
{
    // Define some data points, which should be searched.
    // Each column defines one datapoint.
    Matrix dataPoints(3, 9);
    dataPoints << 1, 2, 3, 1, 2, 3, 1, 2, 3,
                  2, 1, 0, 3, 2, 1, 0, 3, 4,
                  3, 1, 3, 1, 3, 4, 4, 2, 1;

    kdt::KDTreed kdtree(dataPoints);
    kdtree.build();

    // Create a querypoint. We will search for this points nearest neighbors.
    Matrix queryPoints(3, 1);
    queryPoints << 1, 2, 3;

    // Initialize result matrices.
    kdt::KDTreed::Matrix dists; // basically Eigen::MatrixXd
    kdt::KDTreed::MatrixI idx; // basically Eigen::Matrix<Eigen::Index>
    size_t knn = 3;
    kdtree.query(queryPoints, knn, idx, dists);

    // Do something with the results.
    std::cout
        << "Data points:" << std::endl
        << dataPoints << std::endl
        << "Query points:" << std::endl
        << queryPoints << std::endl
        << "Neighbor indices:" << std::endl
        << idx << std::endl
        << "Neighbor distances:" << std::endl
        << dists << std::endl;

    return 0;
}