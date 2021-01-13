/*
 * GPShape.cpp: GPShape functions
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

#include "GPShape.hpp"
#include <algorithm>

// Update the Kxx kernel (3*N, 3*N) with scaled measurements. Add noise along
// the diagonal elements. Mean and variance values must be unscaled to make
// sense
void GPShape::updateKxx(const int& idx){
    // std::cout<<"[GPShape::updateKxx]\n";
    int n = static_cast<int>(gps[idx]->measurements.size());
    if (n == 1)     gps[idx]->Kxx = Emx::Zero(3*n,3*n); // n is the number of training points
    else            gps[idx]->Kxx.conservativeResize(3*n, 3*n); // n is the number of training points

    // 0 to n blocks
    for(int i = 0; i != n; i++) {
        std::vector<double> x1 = {gps[idx]->measurements[i].x, gps[idx]->measurements[i].y}; 
        std::vector<double> x2 = {gps[idx]->measurements[n-1].x, gps[idx]->measurements[n-1].y}; 
        gps[idx]->Kxx.block<3,3>(3*i,3*(n-1)) = kernelBlock(x1, x2); // kernel block
    }

    for(int j = 0; j != n; j++) {
        std::vector<double> x1 = {gps[idx]->measurements[n-1].x, gps[idx]->measurements[n-1].y}; 
        std::vector<double> x2 = {gps[idx]->measurements[j].x, gps[idx]->measurements[j].y}; 
        gps[idx]->Kxx.block<3,3>(3*(n-1),3*j) = kernelBlock(x1, x2); // kernel block
    }
    
    addKernelNoise(gps[idx]->Kxx, gps[idx]->numPrior); // add noise to diag

    if (n==1)   gps[idx]->Kxx_L = getCholeskyL(gps[idx]->Kxx); // update cholesky
    else        updateCholeskyL(idx, gps[idx]->Kxx); // update cholesky
}

void GPShape::addKernelNoise(Emx& kernel, const int& num_prior) {
    // std::cout<<"[GPShape::addKernelNoise]\n";
    int n = kernel.rows()/3; 
    Emx noise = Emx::Zero(3*n,3*n); 
    if ((n -1) < num_prior){
        for (int i = 0; i < 3; i++)
            noise(3*(n -1) + i,3*(n -1) + i) += prior_var[i];
    }
    else {
        for (int i = 0; i < 3; i++)
            noise(3*(n -1) + i,3*(n -1) + i) += var[i];   
    }
    kernel = kernel + noise; // add noise to diag 
}

// Get the Kxstarxstar kernel (3, 3)
void GPShape::getKxStar2(const Emx& x_star, Emx& K_xStar_xStar) const{
    // std::cout<<"[GPShape::getKxStar2]\n";
    int sz = static_cast<int>(x_star.rows()); 
    for (int i = 0; i < sz; i++) {
        std::vector<double> xi = {x_star(i, 0), x_star(i, 1)};
        for (int j = 0; j < sz; j++) {
            std::vector<double> xj = {x_star(j, 0), x_star(j, 1)};
            K_xStar_xStar.block<3,3>(3*i, 3*j) = kernelBlock(xi, xj);
        }
    }
}

// Get the Kxstarxstar kernel (3, 3)
void GPShape::getKxStar2(const std::vector<double>& x, Emx& K_xStar_xStar) const{
    // std::cout<<"[GPShape::getKxStar2]\n";
    K_xStar_xStar = kernelBlock(x, x);
}

// get the Kxstar kernel (3*N, 3) from scaled measurements
void GPShape::getKxStar(const int& idx, const std::vector<double>& x, Emx& KxStar) const{
    // std::cout<<"[GPShape::getKuStar]\n";
    int n = static_cast<int>(gps[idx]->measurements.size());
    KxStar = Emx(3*n, 3);
    for(int i = 0; i != n; i++) {
        std::vector<double> xm = {gps[idx]->measurements[i].x, gps[idx]->measurements[i].y}; 
        KxStar.block<3,3>(3*i, 0) = kernelBlock(xm, x);   // add a 3 x 3 block 
    }
}   

// get the Kxstar kernel (3*N, 3) from scaled measurements
void GPShape::getKxStar(const int& idx, const Emx& x_star, Emx& KxStar) const{
    // std::cout<<"[GPShape::getKuStar]\n";
    // xstar = 3n, 3testlen
    int n = gps[idx]->measurements.size(); // number meas
    KxStar = Emx(3*n, 3);
    for(int i = 0; i != n; i++) {
        std::vector<double> xi = {gps[idx]->measurements[i].x, gps[idx]->measurements[i].y}; 
        for (int j = 0; j < test_side2; j++) {
            std::vector<double> xj = {x_star(j, 0), x_star(j, 1)};
            KxStar.block<3,3>(3*i, 3*j) = kernelBlock(xi, xj);
        }
    }
} 

// get the Kxstar kernel (3*N, 3) from scaled measurements
void GPShape::getKxCheck(const int& idx, const Emx& x_check, Emx& KxCheck) const{
    // std::cout<<"[GPShape::getKxCheck]\n";
    int n = static_cast<int>(gps[idx]->measurements.size()); // number meas
    std::vector<double> testpoint = {x_check(0, 0), x_check(0, 1)};
    for (int i = 0; i < n; i++) {
        std::vector<double> xi = {gps[idx]->measurements[i].x, gps[idx]->measurements[i].y}; 
        KxCheck.block<3,3>(3*i, 0) = kernelBlock(xi, testpoint);
    }
}   

// Add new contact and normal measurements
void GPShape::addMeasurement(const int& idx, const Evx& contact,  const Evx& normal){
    // std::cout<<"[GPShape::addMeasurement]\n";
    // normalize normals
    double n_x = normal.x()/std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2));
    double n_y = normal.y()/std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2)); 
    Contact c(contact[0], contact[1], n_x, n_y);
    gps[idx]->measurements.push_back(c);  
    gps[idx]->numMeasurements++;

    gps[idx]->outputs.conservativeResize(gps[idx]->outputs.rows() + 3, gps[idx]->outputs.cols());
    gps[idx]->outputs(gps[idx]->outputs.rows() - 3, 0) = 0.0;
    gps[idx]->outputs(gps[idx]->outputs.rows() - 2, 0) = gps[idx]->measurements.back().nx;  
    gps[idx]->outputs(gps[idx]->outputs.rows() - 1, 0) = gps[idx]->measurements.back().ny;  

}

bool GPShape::preFilter(const int& idx, const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p) { 
    // std::cout<<"[GPShape::preFilter]\n";
    Evx p_transform = transformPointTo(contact_pt, p); // to local frame
    Evx n_transform = transformPointTo(contact_normal, p);

    double v_pred = 0.0, n_dot_pred = 0.0, d_pred = 0.0;
    // screen new measurement
    if (gps[idx]->numMeasurements > gps[idx]->numPrior - 1) {
        testMeasurement(idx, p_transform, n_transform, v_pred, n_dot_pred, d_pred);
    } else {
        return true;
    }

    return (v_pred > vThresh);
}

bool GPShape::proximityCheck(const int& idx, const Evx& x) const {
    // std::cout<<"[GPShape::proximityCheck]\n";
    double dist2center = sqrt((gps_centers[idx][0] - x(0))*(gps_centers[idx][0] - x(0)) + (gps_centers[idx][1] - x(1))*(gps_centers[idx][1] - x(1)));
    if (dist2center < gp_radius)    return true;
    else                            return false;
}


void GPShape::updateThread(const int& start_id, const int& end_id,  const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p, std::vector<bool>& isUpdated) {
    // std::cout<<"[GPShape::updateThread]\n";
    for (int i = start_id; i < end_id; i++) {
        isUpdated[i]= updateLocal(i, contact_pt, contact_normal, p);
    }
}


bool GPShape::update(const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p) { 
    // std::cout<<"[GPShape::update]\n";
    Timer update_timer; 

    isUpdated = std::vector<bool>(gp_count, false);
    
    if(localExperts) {
        #if 1
            int n_threads = std::thread::hardware_concurrency();
            int n_threads_used = n_threads;
            if (gp_count < n_threads) n_threads_used = test_side2;
            else                    n_threads_used = n_threads;

            std::thread *threads = new std::thread[n_threads_used];

            int n_excess = gp_count % n_threads_used; // excess test points 
            int batch_sz = gp_count / n_threads_used; // size per thread
            int batch_start = 0;

            // allot batch size + 1 for the excess 
            for(int i = 0; i < n_excess; ++i){
                threads[i] = std::thread(&GPShape::updateThread, this, batch_start,
                        batch_start + batch_sz + 1,  contact_pt, contact_normal, p, std::ref(isUpdated));
                batch_start += batch_sz + 1;
            }

            // allot batch size for the rest 
            for (int i = n_excess; i < n_threads_used; ++i){
                threads[i] = std::thread(&GPShape::updateThread, this, batch_start,
                    batch_start + batch_sz,  contact_pt, contact_normal, p, std::ref(isUpdated));
                batch_start += batch_sz;
            }

            // wait to complete 
            for (int i = 0; i < n_threads_used; ++i){
                threads[i].join();
            }

            delete [] threads;
        #else  
            for (int i = 0; i < gp_count; i++) {
                update(i, contact_pt, contact_normal, p);
            }
        #endif
    }
    else {
        isUpdated[0] = updateLocal(0, contact_pt, contact_normal, p);
    }

    bool updated = !std::all_of(isUpdated.begin(), isUpdated.end(), [](bool v) { return !v; });
    if(localExperts)    std::cout << "[GPShape::update] "<< update_timer.elapsed() << " secs (updated: "<<updated<<") LocalExperts\n";
    else                std::cout << "[GPShape::update] "<< update_timer.elapsed() << " secs (updated: "<<updated<<")\n";

    update_time.push_back(update_timer.elapsed());
    return updated; // all false
}

// Run all update steps given new observations
bool GPShape::updateLocal(const int& idx, const Evx& contact_pt,  const Evx& contact_normal,  const Evx& p, const bool& isPrior) { 

    gps[idx]->currPose = p; 
    Evx p_transform = transformPointTo(contact_pt, gps[idx]->currPose); // to local frame
    Evx n_transform = transformPointTo(contact_normal, gps[idx]->currPose);
    
    // not belonging to local GP, and isn't prior data
    if(localExperts && (!isPrior) && (!proximityCheck(idx, p_transform))) {
        return false;
    }

    if(isPrior) {
        // empty GP
        // std::cout<<"[GPShape::update] add Prior "<<idx<<std::endl;
        addMeasurement(idx, p_transform, n_transform); 
        updateKxx(idx); 
        updateAlpha(idx);       
    }
    else if(preFilter(idx, contact_pt, contact_normal, p) || noPreFilter) {
        // useful measurements as inducing pts
        // std::cout<<"[GPShape::update] add Measurement "<<idx<<std::endl;
        addMeasurement(idx, p_transform, n_transform); 
        updateKxx(idx); 
        updateAlpha(idx);
    }
    else {
        return false;
    }

    return true;
}


// Run all update steps given new observations
void GPShape::addPrior(const double& r, const int& n, const Evx& p) { 
    // std::cout<<"[GPShape::addPrior]\n";
    std::vector<Evx> prior_pts, prior_normals; 
    for (int i = 0; i < n; i++) {
        const double rad = ((double)i / (double) n) * M_PI * 2;
        Evx pt(2); pt << cos(rad) * r, sin(rad) * r;
        Evx p_transform = transformPointFrom(pt, p);
        Evx contact_pt(2), contact_normal(2);

        contact_pt <<p_transform[0], p_transform[1] ;
        contact_normal << p_transform[0], p_transform[1] ;

        prior_pts.push_back(contact_pt);
        prior_normals.push_back(contact_normal);
    }

    if(localExperts) {
        for (int j = 0; j < gp_count; j++) {
            for (int i = 0 ; i < prior_pts.size(); i++) {
                // prior added to all local experts
                if(updateLocal(j, prior_pts[i], prior_normals[i], p, true)) {
                    gps[j]->numPrior++;
                }
            }
        }
    } else {
            for (int i = 0 ; i < prior_pts.size(); i++) {
                // prior added to all local experts
                if(updateLocal(0, prior_pts[i], prior_normals[i], p, true)) {
                    gps[0]->numPrior++;
                }
            }
    }
}

void GPShape::updateAlpha(const int& idx) {
    // std::cout<<"[GPShape::updateAlpha] "<<idx<<"\n";
    if (localExperts)   gps[idx]->alpha = solveCholeskyWithL(gps[idx]->Kxx_L, gps[idx]->outputs);
    else        gps[0]->alpha = solveCholesky(gps[0]->Kxx, gps[0]->outputs);
}

Emx GPShape::getCholeskyL(const Emx& A) const {
    return Eigen::LLT<Eigen::MatrixXd>(A).matrixL(); 
}

Emx GPShape::solveCholeskyWithL(const Emx& L, const Emx& b) const {
    Emx L_x = L.triangularView<Eigen::Lower>().solve(b);
    return (L.triangularView<Eigen::Lower>().transpose()).solve(L_x);
}

Emx GPShape::solveCholesky(const Emx& A, const Emx& b) const {
    Emx L = getCholeskyL(A); 
    return solveCholeskyWithL(L, b);
}

void GPShape::updateCholeskyL(const int& idx, const Emx& A) {
    int n = A.rows() - 3;
    Emx A12 = A.block(0,n,n,3);
    Emx A22 = A.block(n, n, 3, 3);

    Emx L12 = gps[idx]->Kxx_L.triangularView<Eigen::Lower>().solve(A12);
    Emx L22 = Eigen::LLT<Eigen::MatrixXd>(A22 - L12.transpose() * L12).matrixL(); 

    gps[idx]->Kxx_L.conservativeResize(gps[idx]->Kxx_L.rows() + 3, gps[idx]->Kxx_L.cols() + 3);
    gps[idx]->Kxx_L.block(n,0,3,n) = L12.transpose(); 
    gps[idx]->Kxx_L.block(n, n, 3, 3) = L22;     
}

// Test meshgrid, get mean and variance matrices, and compute the implicit
// surface. Multithreaded solve.
void GPShape::test(cnt& contour, Emx& m_test, Emx& c_test, Evx& contourVar){
    Timer test_timer; 

    if (M_test.isZero())
        isUpdated = std::vector<bool>(gp_count, true);

    #if 1
        /////// threads 
        int n_threads = std::thread::hardware_concurrency();
        int n_threads_used = n_threads;
        if (test_side2 < n_threads) n_threads_used = test_side2;
        else                    n_threads_used = n_threads;

        std::thread *threads = new std::thread[n_threads_used];

        int n_excess = test_side2 % n_threads_used; // excess test points 
        int batch_sz = test_side2 / n_threads_used; // size per thread
        int batch_start = 0;

        // allot batch size + 1 for the excess 
        for(int i = 0; i < n_excess; ++i){
            threads[i] = std::thread(&GPShape::testThread, this, i, batch_start,
                 batch_start + batch_sz + 1, test_side, std::ref(M_test), std::ref(C_test));
            batch_start += batch_sz + 1;
        }

        // allot batch size for the rest 
        for (int i = n_excess; i < n_threads_used; ++i){
            threads[i] = std::thread(&GPShape::testThread, this, i, batch_start,
                batch_start + batch_sz, test_side, std::ref(M_test), std::ref(C_test));
            batch_start += batch_sz;
        }

        // wait to complete 
        for (int i = 0; i < n_threads_used; ++i){
            threads[i].join();
        }

        delete [] threads;
    #elif 0
        M_test  = Emx::Zero (test_side, test_side); 
        C_test  = Emx::Zero (test_side, test_side); 
        testThread(-1, 0, test_side2, test_side, M_test, C_test); 
    #else
        int n = static_cast<int>(inducing_points.size()); // n is the number of training points 
        Emx x_star = Emx::Zero(test_side2, 2);

        for (int i = 0; i < test_side2; ++i) {
            int r = i/test_side; int c = i%test_side;
            x_star(i, 0) = gpis_->x_test(r, c); 
            x_star(i, 1) = gpis_->y_test(r, c); 
        }

        M_test  = Emx::Zero (test_side, test_side); 
        C_test  = Emx::Zero (test_side, test_side); 

        Emx K_uStar = Emx::Zero(3*n, 3*test_side2); // n is the number of training points
        Emx K_xStar_xStar = Emx::Zero(3*test_side2,3*test_side2); // n is the number of training points

        getKxStar2(x_star, K_xStar_xStar);
        getKuStar(x_star, K_uStar);  
        Emx M, V, C;

        M =  K_uStar.transpose() * alpha;   // mean equation
        C = K_xStar_xStar -  K_uStar.transpose() * beta;

        for (int i = 0; i < test_side2; i++) {
            int r = i/test_side; int c = i%test_side;
            M_test(r, c) = M(3*i, 0);
            C_test(r, c) = C(3*i, 3*i);
        }
    #endif

    // get GPIS contour
    std::unique_ptr<CContourMap> map = std::make_unique<CContourMap>();
    map->generate_levels(0, 0, 1);
    map->contour(M_test, xx, yy); 
    map->consolidate(); 
    map->getBestShapeContour(contour, testLim); // single thresholded contour, testLim*scale = half side
    contourVar = getContourVariance(contour, C_test);
    std::cout << "[GPShape::test] "<< test_timer.elapsed() << " secs"<<"\n";
    test_time.push_back(test_timer.elapsed());

    m_test = M_test; c_test = C_test;
}

Evx GPShape::getContourVariance(const cnt& contour, const Emx& c_test) const {
    int test_sz = static_cast<int>(x_test.size1()); // number of test elements 
    int contour_sz = static_cast<int>(contour.size());
    Evx fVar(contour_sz);
    for (int i = 0; i < contour_sz; ++i) {
        std::vector<double> p = {contour[i][0], contour[i][1]}; // scale back
        std::vector<int> p_grid = {(int)(p[0]/res) + test_sz/2,(int)(p[1]/res) + test_sz/2}; // scale back
        fVar(i) = C_test(p_grid[0], p_grid[1]);
    }
    return fVar;
}

// Threads take subset of entire test points
void GPShape::testThread(const int& thread_id, const int& start_id,
                     const int& end_id, const int& sz, Emx& m_, Emx& c_) const {
    // printf ("testThread: %d start_id: %d end_id: %d sz: %d\n", thread_id, start_id, end_id, sz);

    if (localExperts) {
        Emx K_xStar; //= Emx::Zero(3*n, 3); // n is the number of training points
        Emx K_xStar_xStar = Emx::Zero(3,3); // n is the number of training points

        std::vector<double> x_star;
        Emx M, V, C;

        for (int i = start_id; i < end_id; ++i) {
            int r = i/sz; int c = i%sz;
            x_star = {x_test(r, c), y_test(r, c)};
            getKxStar2(x_star, K_xStar_xStar);
            
            double m__ = 0.0; double c__ = 0.0; 
            for (int j = 0; j < nnGPValid[i].size(); ++j) {
                int active_gp_id = nnGPValid[i][j]; 

                if (!isUpdated[active_gp_id]) {
                    m__ += m_(r,c);
                    c__ += c_(r,c);
                    continue;
                }
                    
                getKxStar(active_gp_id, x_star, K_xStar);     
                M =  K_xStar.transpose() *  gps[active_gp_id]->alpha;   // mean equation
                m__ += M(0, 0);

                V = gps[active_gp_id]->Kxx_L.triangularView<Eigen::Lower>().solve(K_xStar);
                C = K_xStar_xStar - V.transpose()*V;
                c__ += C(0, 0);
            }
            m_(r,c) = m__/nnGPValid[i].size();
            c_(r,c) = c__/nnGPValid[i].size();  // scale variance   
        }
    }
    else {
        Emx K_xStar = Emx::Zero(3*gps[0]->measurements.size(), 3); // n is the number of training points
        Emx K_xStar_xStar = Emx::Zero(3,3); // n is the number of training points

        std::vector<double> x_star;
        Emx M, V, C;

        for (int i = start_id; i < end_id; ++i) {
            int r = i/sz; int c = i%sz;

            x_star = {x_test(r, c), y_test(r, c)};

            getKxStar2(x_star, K_xStar_xStar);
            getKxStar(0, x_star, K_xStar);    

            M =  K_xStar.transpose() * gps[0]->alpha;   // mean equation
            m_(r,c) = M(0, 0);
            
            V = gps[0]->Kxx_L.triangularView<Eigen::Lower>().solve(K_xStar);
            C = K_xStar_xStar - V.transpose()*V;
            c_(r,c) = C(0, 0);  // scale variance    
        }
    }
}

// Screen measurement and get expected variance, expected-measured normal dot
// product, expected d
void GPShape::testMeasurement(const int& idx, const Evx& p, const Evx& n, 
                                double& v, double& n_dot, double& d) const{

    if (localExperts) {
        Emx x_t(1, 2); x_t(0,0) = p(0); x_t(0,1) = p(1);
        Emx K_xTest = Emx::Zero(3*gps[idx]->measurements.size(), 3); // n is the number of training points
        Emx K_xTest_xTest = Emx::Zero(3,3); // n is the number of training points

        getKxStar2(x_t, K_xTest_xTest);
        getKxCheck(idx, x_t, K_xTest); 

        Emx M =  K_xTest.transpose() * gps[idx]->alpha;   // mean equation
        d = std::abs(M(0, 0)); // scaled sdf value
        Evx n_pred(2); n_pred << M(1, 0)/std::sqrt(std::pow(M(1, 0), 2) + std::pow(M(2, 0), 2)),
                                M(2, 0)/std::sqrt(std::pow(M(1, 0), 2) + std::pow(M(2, 0), 2));
        Evx n_meas(2); n_meas << n(0)/std::sqrt(std::pow(n(0), 2) + std::pow(n(1), 2)),
                                n(1)/std::sqrt(std::pow(n(0), 2) + std::pow(n(1), 2));
        n_dot = n_pred.dot(n_meas); // normal dot product

        Emx V = gps[idx]->Kxx_L.triangularView<Eigen::Lower>().solve(K_xTest);
        Emx C = K_xTest_xTest - V.transpose()*V;
        v = C(0, 0);  // scale variance    
    }
    else {
        Emx x_t(1, 2); x_t(0,0) = p(0); x_t(0,1) = p(1);
        Emx K_xTest = Emx::Zero(3*gps[0]->measurements.size(), 3); // n is the number of training points
        Emx K_xTest_xTest = Emx::Zero(3,3); // n is the number of training points

        getKxStar2(x_t, K_xTest_xTest);
        getKxCheck(0, x_t, K_xTest); 

        Emx M =  K_xTest.transpose() * gps[0]->alpha;   // mean equation
        d = std::abs(M(0, 0)); // scaled sdf value
        Evx n_pred(2); n_pred << M(1, 0)/std::sqrt(std::pow(M(1, 0), 2) + std::pow(M(2, 0), 2)),
                                M(2, 0)/std::sqrt(std::pow(M(1, 0), 2) + std::pow(M(2, 0), 2));
        Evx n_meas(2); n_meas << n(0)/std::sqrt(std::pow(n(0), 2) + std::pow(n(1), 2)),
                                n(1)/std::sqrt(std::pow(n(0), 2) + std::pow(n(1), 2));
        n_dot = n_pred.dot(n_meas); // normal dot product

        Emx V = gps[0]->Kxx_L.triangularView<Eigen::Lower>().solve(K_xTest);
        Emx C = K_xTest_xTest - V.transpose()*V;
        v = C(0, 0);  // scale variance    
    }
}

// Generate the (3, 3) kernel block using specified kernel equations
// incorporating gradient information from normals. For (x1, x1) we only update
// the diagonal terms
Emx GPShape::kernelBlock(const std::vector<double>& x1, const std::vector<double>& x2) const{
    if (kernel == "thinplate")      return kernelBlockThinPlate(x1, x2);
    else if (kernel == "gaussian")  return kernelBlockGaussian(x1, x2);
    else if (kernel == "matern")    return kernelBlockMatern(x1, x2);
    else                            return kernelBlockThinPlate(x1, x2);
}

// Use thin plate kernel
Emx GPShape::kernelBlockThinPlate(const std::vector<double>& x1, 
                                  const std::vector<double>& x2) const{
    Emx block =  Emx::Zero(3,3);
    double xyNorm = Norm(x1, x2);
    if (xyNorm > 1e-6) {
        // x1 != x2
        block(0,0) = 2.0 * std::pow(xyNorm, 3) - 3.0 * R * std::pow(xyNorm, 2) + std::pow(R, 3); //  cov(di, dj) 
        block(0,1) = 2.0 * 3.0 * xyNorm * (-x1[0] + x2[0]) - 3.0 * R * 2.0 * (-x1[0] + x2[0]); // cov(di, w1j)
        block(0,2) = 2.0 * 3.0 * xyNorm * (-x1[1] + x2[1]) - 3.0 * R * 2.0 * (-x1[1] + x2[1]); // cov(di, w2j)

        block(1,0) = -block(0,1); // cov(w1i, dj)
        block(1,1) = 2.0 * 3.0 * ((-x1[0] + x2[0])/xyNorm) * (x1[0] - x2[0]) + 2.0 * 3.0 * xyNorm * (-1.0) - 3.0 * R * 2.0 * (-1.0); // cov(w1i, w1j)
        block(1,2) = 2.0 * 3.0 * ((-x1[1] + x2[1])/xyNorm) * (x1[0] - x2[0]); // cov(w1i, w2j)

        block(2,0) = -block(0,2); // cov(w2i, dj)
        block(2,1) =  block(1,2); // cov(w2i, w1j)
        block(2,2) = 2.0 * 3.0 * ((-x1[1] + x2[1])/xyNorm) * (x1[1] - x2[1]) + 2.0 * 3.0 * xyNorm * (-1.0) - 3.0 * R * 2.0 * (-1.0); // cov(w2i, w2j)
    } else {
        // x1 == x2 diag
        block(0,0) = 2.0 * std::pow(xyNorm, 3) - 3.0 * R * std::pow(xyNorm, 2) + std::pow(R, 3); 
        block(1,1) = 6.0 * R; 
        block(2,2) = 6.0 * R; 
    }
    return block;
}

// Use squared exponential kernel
Emx GPShape::kernelBlockGaussian(const std::vector<double>& x1, 
                                 const std::vector<double>& x2) const{
    Emx block =  Emx::Zero(3,3);
    double l2_inv = 1.0 / std::pow(l, 2);

    double diffX = x1[0] - x2[0];
    double diffY = x1[1] - x2[1];

    block(0,0) = 1.0;

    block(1,1) = l2_inv * (1 - std::pow(diffX,2) * l2_inv);
    block(2,2) = l2_inv * (1 - std::pow(diffY,2) * l2_inv);
    block(1,2) = -std::pow(l2_inv,2) * diffX * diffY;

    block(0,1) = -l2_inv * diffX;
    block(0,2) = -l2_inv * diffY;

    block(1,0) = block(0,1);
    block(2,0) = block(0,2);
    block(2,1) = block(1,2);

    block *= std::pow(sigma, 2) * std::exp(-0.5 * Norm(x1, x2) * l2_inv);

    return block;
}

// Use matern kernel
Emx GPShape::kernelBlockMatern(const std::vector<double>& x1, 
                               const std::vector<double>& x2) const{
    Emx block =  Emx::Zero(3,3);
    double d = Norm(x1, x2);
    double root5 = std::sqrt(5);
    double rho2 = std::pow(rho,2);

    double diffX = x1[0] - x2[0];
    double diffY = x1[1] - x2[1];

    block(0,0) = 1.0 + (root5 * d / rho) + 4.0 * std::pow(d,2) / (3.0 * rho2);
    block(1,1) = -5.0 / (3.0 * std::pow(rho,4)) * 
        (5.0 * std::pow(diffX,2) - rho2 - rho * root5 * d);
    block(2,2) = -5.0 / (3.0 * std::pow(rho,4)) * 
        (5.0 * std::pow(diffY,2) - rho2 - rho * root5 * d);

    double temp = -5.0 * (root5 * d + rho) / (3.0 * std::pow(rho,3));
    block(0,1) = diffX * temp;
    block(0,2) = diffY * temp;
    block(1,2) = -25.0 / (3.0 * std::pow(rho,4)) * diffX * diffY;

    block(1,0) = block(0,1);
    block(2,0) = block(0,2);
    block(2,1) = block(1,2);

    block *= std::pow(sigma, 2) * std::exp(-root5 * d / rho);

    return block;
}

void GPShape::buildAndQueryTree() {
    Eigen::MatrixXd dataPoints(2, gp_count);
    for (int i = 0; i < gp_count; i++)
        dataPoints.col(i) = Eigen::VectorXd::Map(&gps_centers[i][0],gps_centers[i].size());
            
    kdt::KDTreed kdtree(dataPoints);
    kdtree.setMaxDistance(gp_radius);
    kdtree.build();

    Eigen::MatrixXd queryPoints(2, test_side2);

    for (int i = 0; i < test_side2; i++) {
        int r = i/test_side; int c = i%test_side;
        queryPoints(0, i) = x_test(r, c);
        queryPoints(1, i) = y_test(r, c);
    }

    // Initialize result matrices.
    kdt::KDTreed::Matrix dists; // basically Eigen::MatrixXd
    kdt::KDTreed::MatrixI idx; // basically Eigen::Matrix<Eigen::Index>
    size_t knn = 3; // top 3 nearest neighbours 
    kdtree.query(queryPoints, knn, nnGP, dists);

    // copy valid GPs
    nnGPValid = std::vector <std::vector <int> >(nnGP.cols(), vector<int>());
    for (int i = 0; i < nnGP.cols(); i++) {
        for (int j = 0; j < nnGP.rows(); j++) {
            int gp_idx = (int) nnGP(j, i);
            if (gp_idx == -1)   break;
            else                nnGPValid[i].push_back(gp_idx);
        }
    }
}

void GPShape::initGPs() {
    int length = (int)sqrt(gp_count); 
    double dist = 2*testLim/((double)length + 1.0);
    gp_radius = 2*dist; // 2 times center spacing 
    for (int i = 0; i < gp_count; i++) {
        gps.push_back(new GP(i, gp_count, testLim));
        gps_centers.push_back(gps[i]->getGPCenter());
    }
}

// Transform contour to pose
void GPShape::transformContour(cnts& contours, const Evx& p) const{
    for (uint i = 0; i < contours.size(); ++i) {
        for(uint j = 0; j<contours[i].size(); j++) {
            Evx pt(2); pt << contours[i][j][0], contours[i][j][1];
            Evx p_transform_back = transformPointFrom(pt, p);
            contours[i][j][0] = p_transform_back[0];
            contours[i][j][1] = p_transform_back[1];
        } 
    }
}

// From local to global frame
Evx GPShape::transformPointFrom(const Evx& pt, const Evx& pose2) const{
    Evx t_pt(2); t_pt << cos(pose2[2]) * pt[0] - sin(pose2[2]) * pt[1], sin(pose2[2]) * pt[0] + cos(pose2[2]) * pt[1]; 
    t_pt = t_pt + pose2.head(2); 
    return t_pt;
}

// From global to local frame
Evx GPShape::transformPointTo(const Evx& pt, const Evx& pose2) const{
    Evx d = pt - pose2.head(2); 
    Evx t_pt(2); t_pt << cos(pose2[2]) * d[0] + sin(pose2[2]) * d[1], -sin(pose2[2]) * d[0] + cos(pose2[2]) * d[1]; 
    return t_pt;
}

// Euclidean distance norm Norm(x - y)
double GPShape::Norm( const std::vector<double>& x,  const std::vector<double>& y) const {
    //calculating number to square in next step
	double a_ = x[0] - y[0]; 
	double b_ = x[1] - y[1];
    //calculating Euclidean distance
	return std::sqrt(std::pow(a_, 2) + std::pow(b_, 2)); 
}

