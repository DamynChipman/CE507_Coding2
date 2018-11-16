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

// ----- Project Files -----
#include "BSpline.h"
#include "DBSpline.h"

using namespace Eigen;

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> tri;

Eigen::VectorXf FE1D(int* ID, int** IEN, int** LM, int p, int Ne, int NINT) {
    
    // Initialize matrices
    std::vector<tri> coefs;
    Eigen::VectorXf F(Ne);
    for (int i = 0; i < Ne; i++) {
        F(i) = 0;
    }
    
    // Necessary Variables
    
    
    // LOOP GOES HERE
    for (int e = 0; e < Ne; e++) {
        for (int i = 0; i < NINT; i++) {
            // Compute Bpi
            BSpline Bpi(p, NINT);
            Eigen::VectorXf Bpi_vec = Bpi.getVector();
            
            // Compute Np
            Eigen::VectorXf Np_vec = C_P(p, e)*Bpi_vec;
            
            // Compute DBpi
            DBSpline DBpi(p, NINT);
            Eigen::VectorXf DBpi_vec = DBpi.getVector();
            
            // Compute DNp
            Eigen::VectorXf DNp_vec = C_P(p, e)*DBpi_vec;
            
            // Compute x(xi_i)
            
            // Compute f_i = f(x(xi_i))
            
            for (int a = 0; a < p+1; a++) {
                for (int b = 0; b < p+1; b++) {
                    // Compute and update k_e_ab
                    
                }
                // Compute and update f_e_a
                
            }
        }
        for (int a = 0; a < p+1; a++) {
            for (int b = 0; b < p+1; b++) {
                if (LM[a][e] > -1) {
                    // Compute and update K_LM[a][e]_LM[a][e]
                    
                }
            }
            if (LM[a][e] > -1) {
                // Compute and update F_LM[a][e]
                
            }
        }
    }
    
    // --------------
    
    SpMat K(Ne,Ne);
    K.setFromTriplets(coefs.begin(), coefs.end());
    Eigen::SimplicialCholesky<SpMat> chol(K);
    Eigen::VectorXf d = chol.solve(F);
    return d;
}



#endif /* FE1D_h */
