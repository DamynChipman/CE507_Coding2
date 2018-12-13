//   Created by Damyn Chipman on 10/26/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding2

// ----- GLOBAL VARIABLES -----
bool VERBOSE = true;                 // Option for printing information
bool VERIFY = false;                 // Option for verifying cases
int NNODES = 20;                     // Number of nodes per elements
int NELEM;                           // Number of elements
int N[4] = {10, 100, 1000, 10000};   // Number of elements array

// ----- Necessary Files -----
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// ----- Helper Function Declarations -----
Eigen::Matrix3f C_P2(int e);
Eigen::Matrix4f C_P3(int e);
Eigen::MatrixXf C_P(int p, int e);
float f_x(float x);
float u_actual(float x);

// ----- Project Files -----
#include "Verifications.h"
#include "Domain1D.h"
#include "BSpline.h"
#include "DBSpline.h"
#include "FE1D.h"

// ----- Namespaces -----
using namespace std;
using namespace Eigen;

// ----- MAIN FUNCTION -----
/**
 * @function main
 * @brief Main driver for CE507_Coding2
 * @discussion Implements a BSpline basis function with extraction operators to Bezier curves.
 * @returns 1 : int : Program ran successfully
 * @returns 9 : int : Program terminated prematurely (user through debugging/errors)
 */
int main() {
    
    int n = 0;             // Number of elements index
    int p = 2;             // BSpline order
    
    // Verification Testing ========================================================
    
    if (VERIFY) { cout << "   running verifications..." << endl; RunVerifications(); }
    
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
    for (int e = 0; e < NE; e++) {
        
    }
//    float u_FE[NE];
//    for (int i = 0; i < NE; i++) {
//        u_FE[i] = 0.0;
//        for (int a = 0; a < p+1; a++) {
//            u_FE[i] += d(i+a);
//        }
//    }
    
    
    return 0;
}

// ----- Helper Function Definitions -----

Eigen::Matrix2f C_P1(int e) {
    Eigen::Matrix2f C;
    
    C(0,0) = 1.0; C(0,1) = 1.0;
    C(1,0) = 1.0; C(1,1) = 1.0;
    
    return C;
}

/**
 * @function C_P2
 * @brief Extraction operator matrix, P = 2
 * @param e : int : Element number
 * @returns C : Matrix3f : Extraction matrix
 */
Eigen::Matrix3f C_P2(int e) {
    Eigen::Matrix3f C;
    
    if (e == 0) {
        C(0,0) = 1.0;   C(0,1) = 0.0;   C(0,2) = 0.0;
        C(1,0) = 0.0;   C(1,1) = 1.0;   C(1,2) = 0.5;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 0.5;
    }
    else if (e == NELEM-1) {
        C(0,0) = 0.5;   C(0,1) = 0.0;   C(0,2) = 0.0;
        C(1,0) = 0.5;   C(1,1) = 1.0;   C(1,2) = 0.0;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 1.0;
    }
    else {
        C(0,0) = 0.5;   C(0,1) = 0.0;   C(0,2) = 0.0;
        C(1,0) = 0.5;   C(1,1) = 1.0;   C(1,2) = 0.5;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 0.5;
    }
    return C;
}

/**
 * @function C_P3
 * @brief Extraction operator matrix, P = 3
 * @param e : int : Element number
 * @returns C : Matrix4f : Extraction matrix
 */
Eigen::Matrix4f C_P3(int e) {
    Eigen::Matrix4f C;
    
    if (e == 0) {
        C(0,0) = 1.0;   C(0,1) = 0.0;   C(0,2) = 0.0;   C(0,3) = 0.0;
        C(1,0) = 0.0;   C(1,1) = 1.0;   C(1,2) = 0.5;   C(1,3) = 0.25;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 0.5;   C(2,3) = 7.0/12.0;
        C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = 1.0/6.0;
    }
    else if (e == 1) {
        C(0,0) = 0.25;     C(0,1) = 0.0;     C(0,2) = 0.0;     C(0,3) = 0.0;
        C(1,0) = 7.0/12.0; C(1,1) = 2.0/3.0; C(1,2) = 1.0/3.0; C(1,3) = 1.0/6.0;
        C(2,0) = 1.0/6.0;  C(2,1) = 1.0/3.0; C(2,2) = 2.0/3.0; C(2,3) = 2.0/3.0;
        C(3,0) = 0.0;      C(3,1) = 0.0;     C(3,2) = 0.0;     C(3,3) = 1.0/6.0;
    }
    else if (e == NELEM-2) {
        C(0,0) = 1.0/6.0; C(0,1) = 0.0;     C(0,2) = 0.0;     C(0,3) = 0.0;
        C(1,0) = 2.0/3.0; C(1,1) = 2.0/3.0; C(1,2) = 1.0/3.0; C(1,3) = 1.0/6.0;
        C(2,0) = 1.0/6.0; C(2,1) = 1.0/3.0; C(2,2) = 2.0/3.0; C(2,3) = 7.0/12.0;
        C(3,0) = 0.0;     C(3,1) = 0.0;     C(3,2) = 0.0;     C(3,3) = 1.0/4.0;
    }
    else if (e == NELEM-1) {
        C(0,0) = 1.0/6.0;  C(0,1) = 0.0;   C(0,2) = 0.0;   C(0,3) = 0.0;
        C(1,0) = 7.0/12.0; C(1,1) = 0.5;   C(1,2) = 0.0;   C(1,3) = 0.0;
        C(2,0) = 1.0/4.0;  C(2,1) = 0.5;   C(2,2) = 1.0;   C(2,3) = 0.0;
        C(3,0) = 0.0;      C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = 1.0;
    }
    else {
        C(0,0) = 1.0/6.0; C(0,1) = 0.0;     C(0,2) = 0.0;     C(0,3) = 0.0;
        C(1,0) = 2.0/3.0; C(1,1) = 2.0/3.0; C(1,2) = 1.0/3.0; C(1,3) = 1.0/6.0;
        C(2,0) = 1.0/6.0; C(2,1) = 1.0/3.0; C(2,2) = 2.0/3.0; C(2,3) = 2.0/3.0;
        C(3,0) = 0.0;     C(3,1) = 0.0;     C(3,2) = 0.0;     C(3,3) = 1.0/6.0;
    }
    return C;
}

/**
 * @function C_P
 * @brief Wrapper function for extraction matrices
 * @param p : int : Order
 * @param e : int : Element number
 * @returns C : MatrixXf : Extraction matrix
 */
Eigen::MatrixXf C_P(int p, int e) {
    if (p == 2) { return C_P2(e); }
    else if (p == 1) { return C_P1(e); }
    else { return C_P3(e); }
}

/**
 * @function f_x
 * @brief Forcing function
 * @param x : float : Independent variable
 * @returns x^2 : float
 */
float f_x(float x) { return x*x; }
