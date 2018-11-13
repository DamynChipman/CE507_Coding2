//   Created by Damyn Chipman on 10/31/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   FE1D.h
//   PROJECT:   CE507_Coding2

#ifndef FE1D_h
#define FE1D_h

#include <vector>
#include <Eigen>

using namespace Eigen;

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> T;

Eigen::VectorXf FE1D(int p, int Ne, int Nint) {
    
    // Initialize matrices
    std::vector<T> coefs;
    Eigen::VectorXf F(Ne);
    for (int i = 0; i < Ne; i++) {
        F(i) = 0;
    }
    
    // LOOP GOES HERE
    for (int e = 0; e++; e < Ne) {
        for (int i = 0; i++; i < Nint) {
            
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
