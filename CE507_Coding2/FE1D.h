//   Created by Damyn Chipman on 10/31/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   FE1D.h
//   PROJECT:   CE507_Coding2

#ifndef FE1D_h
#define FE1D_h

// ----- Necessary Files -----
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// ----- Project Files -----
#include "BSpline.h"
#include "DBSpline.h"

using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> triplet;

//Eigen::VectorXf FE1D(int** LM, int p, float (*k_e_ab)(int,int), float (*f_e_a)(int,int), int NE, int NINT, Eigen::VectorXf w, float del_e) {
//Eigen::VectorXf FE1D(int** LM, int p, int NE, int NINT, float (*f)(float), float (*k_update)(int,int,int,int), float (*f_update)(int,int,int), Eigen::VectorXf knotVector, Eigen::VectorXf w, float del_e) {
Eigen::VectorXf FE1D(int** LM, int p, int NE, int NINT, float (*f)(float), Eigen::VectorXf knotVector, int NKNOTS, Eigen::VectorXf w, Eigen::VectorXf intPoints, float del_e) {

    // Initialize matrices
    std::vector<triplet> coefs;
    Eigen::VectorXf F(NE);
    for (int i = 0; i < NE; i++) {
        F(i) = 0;
    }
    
    // Necessary Variables
    Eigen::MatrixXf k_ab(p+1,p+1);
    Eigen::VectorXf f_a(p+1);
    
    // ----- Begin Assembly Algorithm -----
    cout << "---------- BEGINNING ASSEMBLY ALGORITHM ----------" << endl;
    for (int e = 0; e < NE; e++) {
        cout << "  Element: " << e << endl << endl;
        
        for (int i = 0; i < NINT; i++) {
            cout << "    Integration point index: " << i << endl << endl;
            // Compute Bpi
            BSpline Bpi(p, intPoints(i));
            Eigen::VectorXf Bpi_vec = Bpi.getVector();
            cout << "      B_array: [";
            for(int n = 0; n < NINT; n++) { cout << Bpi_vec(n) << " ";}
            cout << "]" << endl;
            
            // Compute Np
            Eigen::VectorXf Np_vec = C_P(p, e)*Bpi_vec;
            cout << "      N_array: [";
            for(int n = 0; n < NINT; n++) { cout << Np_vec(n) << " ";}
            cout << "]" << endl;
            
            // Compute DBpi
            DBSpline DBpi(p, intPoints(i));
            Eigen::VectorXf DBpi_vec = DBpi.getVector();
            cout << "      dB_array: [";
            for(int n = 0; n < NINT; n++) { cout << DBpi_vec(n) << " ";}
            cout << "]" << endl;
            
            // Compute DNp
            Eigen::VectorXf DNp_vec = C_P(p, e)*DBpi_vec;
            cout << "      dN_array: [";
            for(int n = 0; n < NINT; n++) { cout << DNp_vec(n) << " ";}
            cout << "]" << endl;
            
            // Compute x(xi_i)
            float x_xi = 0.0;
            for (int A = 0; A < NKNOTS - p - 1; A++) {
                
                // Calc x_A
                float x_A = 0.0;
                for (int i = 0; i < A + p; i++) {
                    x_A = x_A + (1.0/p)*knotVector(A+i);
                }
                
                // Calc x
                x_xi = x_xi + x_A*Np_vec(A);
            }
            cout << "      x: " << x_xi << endl;
            
            // Compute f_i = f(x(xi_i))
            float f_i = f(x_xi);
            cout << "      f(x): " << f_i << endl;
            
            for (int a = 0; a < p+1; a++) {
                for (int b = 0; b < p+1; b++) {
                    // Compute and update k_e_ab
                    k_ab(a,b) = k_ab(a,b) + DNp_vec(a)*DNp_vec(b)*(2.0/del_e)*w(i);
                    //k_ab(a,b) = k_ab(a,b) + k_update(e,i,a,b);
                }
                // Compute and update f_e_a
                f_a(a) = f_a(a) + Np_vec(a)*f_i*(del_e/2)*w(i);
                //f_a(a) = f_a(a) + f_update(e,i,a);
            }
        }
        for (int a = 0; a < p+1; a++) {
            for (int b = 0; b < p+1; b++) {
                if (LM[a][e] > -1) {
                    // Compute and update K_LM[a][e]_LM[a][e]
                    coefs.push_back(triplet(LM[a][e],LM[a][e],k_ab(a,b)));
                }
            }
            if (LM[a][e] > -1) {
                // Compute and update F_LM[a][e]
                F(LM[a][e]) = F(LM[a][e]) + f_a(a);
            }
        }
    }
    // ----- End Assembly Algorithm -----
    
    SpMat K(NE,NE);
    K.setFromTriplets(coefs.begin(), coefs.end());
    Eigen::SimplicialCholesky<SpMat> chol(K);
    Eigen::VectorXf d = chol.solve(F);
    return d;
}

#endif /* FE1D_h */
