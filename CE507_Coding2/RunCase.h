//   Created by Damyn Chipman on 12/13/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   RunCase.h
//   PROJECT:   CE507_Coding2

#ifndef RunCase_h
#define RunCase_h

#include "FE1D.h"

using namespace std;

float RunCase(int n, int p) {
    // ID, IEN, and LM Arrays ======================================================
    cout << "   generating array mappings..." << endl;
    
    // ID array
    int* ID = new int[N[n]];
    for (int i = 0; i < N[n]+1; i++) {
        ID[i] = i;
    }
    ID[N[n]] = -1;
    
    // IEN array
    int** IEN = new int*[p+1];
    for (int i = 0; i < p+1; i++) {
        IEN[i] = new int[N[n]];
    }
    for (int i = 0; i < p+1; i++) {
        for (int j = 0; j < N[n]; j++) {
            if (i == 0) { IEN[i][j] = ID[j]; }
            else { IEN[i][j] = ID[j]+1; }
        }
    }
    
    // LM array
    int** LM = new int*[p+1];
    for (int i = 0; i < p+1; i++) {
        LM[i] = new int[N[n]];
    }
    for (int i = 0; i < p+1; i++) {
        for (int j = 0; j < N[n]; j++) {
            LM[i][j] = ID[IEN[i][j]];
            //cout << LM[i][j] << " ";
        }
        //cout << endl;
    }
    if (VERBOSE) {
        cout << "LM array: " << endl;
        for (int i = 0; i < p+1; i++) {
            for (int j = 0; j < 5; j++) {
                cout << LM[i][j] << " ";
            }
            cout << "... ";
            for (int j = N[n]-5; j < N[n]; j++) {
                cout << LM[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    // Run Assembly Algorithm FE1D ====================================================
    cout << "   calculating coefs..." << endl;
    int NE = N[n];                                        // Number of elements
    NELEM = NE;                                           // Global scope
    float leftBound = 0.0;                                // Left domain bound
    float rightBound = 1.0;                               // Right domain bound
    float del_e = (rightBound - leftBound)/(NE);          // Element spacing
    int NINT = 3;                                         // Number of integration points for quadrature
    Eigen::VectorXf intPoints(NINT);                      // Vector of integration points
    intPoints(0) = -sqrt(3.0/5.0);                        // Integration point 1
    intPoints(1) = 0.0;                                   // Integration point 2
    intPoints(2) = sqrt(3.0/5.0);                         // Integration point 3
    Eigen::VectorXf w(NINT);                              // Vector of weights for quadrature
    w(0) = 5.0/9.0;                                       // Weight 1
    w(1) = 8.0/9.0;                                       // Weight 2
    w(2) = 5.0/9.0;                                       // Weight 3
    int NKNOTS = (p+1)+(p+1)+(NE)-1;                      // Number of knots
    Eigen::VectorXf knotVector(NKNOTS);                   // Knot vector
    int NDEL_E = 1;                                       // Spacing iterator
    for (int i = 0; i < NKNOTS; i++) {
        if (i < (p+1)) { knotVector(i) = leftBound; }
        if (i > ((p)+NE)) { knotVector(i) = rightBound; }
        if (i > p && i < (p+1+NE)) { knotVector(i) = NDEL_E*del_e; NDEL_E++; }
    }
    if (VERBOSE) {
        cout << "Knot Vector: ";
        for (int i = 0; i < knotVector.size(); i++) {
            cout << knotVector(i) << " ";
        }
        cout << endl;
    }
    
    //Eigen::VectorXf d = FE1D(LM, p, NE, NINT, f_x, k_update, f_update, knotVector, w, del_e);
    Eigen::VectorXf d = FE1D(LM, p, NE, NINT, f_x, knotVector, NKNOTS, w, intPoints, del_e); // Run FE1D
    
    // Generate solutions ===============================================================
    cout << "   generating solutions..." << endl;
    
    
    
    // Calculate error ==================================================================
    cout << "   calculating error..." << endl;
    float error = 0.0;
    for (int e = 0; e < NE; e++) {
        float u_FE[NE];
        float u_ACT = 0.0;
        for (int i = 0; i < NINT; i++) {
            // Compute Bpi
            BSpline Bpi(p, intPoints(i));
            Eigen::VectorXf Bpi_vec = Bpi.getVector();
            
            // Compute Np
            Eigen::VectorXf Np_vec = C_P(p, e)*Bpi_vec;
            
            // Calculate FE solution in element
            u_FE[i] = 0.0;
            for (int a = 0; a < p+1; a++) {
                u_FE[i] += d(i+a)*Np_vec(a);
            }
            
            // Calculate actual solution in element
            // Compute x(xi_i)
            float x_xi = 0.0;
            for (int A = 0; A < p+1; A++) {
                
                // Calc x_A
                float x_A = 0.0;
                for (int j = A+1; j < A + p + 1; j++) {
                    //cout << "S(i) @ " << j << ": " << knotVector(j+e) << endl;
                    x_A = x_A + (1.0/p)*knotVector(j+e);
                }
                //cout << "x_A: " << x_A << endl;
                
                // Calc x
                x_xi = x_xi + x_A*Np_vec(A);
            }
            u_ACT = u_actual(x_xi);
            
            // Sum the error
            error += pow(abs(u_ACT - u_FE[i]),2) * (del_e/2) * w(i);
        }
    }
    VERBOSE = true;
    if (VERBOSE) { cout << "Calculated Error: " << error << endl; }
    return error;
}

#endif /* RunCase_h */
