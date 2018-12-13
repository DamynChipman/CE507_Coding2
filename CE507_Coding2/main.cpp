//   Created by Damyn Chipman on 10/26/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding2

// ----- GLOBAL VARIABLES -----
bool VERBOSE = false;                 // Option for printing information
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
#include "RunCase.h"
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
    
    // Verification Testing ========================================================
    if (VERIFY) { cout << "   running verifications..." << endl; RunVerifications(); }
    
    // Run case
    int n = 3;             // Number of elements index
    int p = 3;             // BSpline order
    float error = RunCase(n,p);
    
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

float u_actual(float x) {
    return (pow(x,2)/2)*(1 - pow(x,2));
}
