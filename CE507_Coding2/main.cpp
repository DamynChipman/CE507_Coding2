//   Created by Damyn Chipman on 10/26/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding2

// ----- Necessary Files -----
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Domain1D.h"
#include "BSpline.h"
#include "FE1D.h"

// ----- Namespaces -----
using namespace std;
using namespace Eigen;

// ----- GLOBAL VARIABLES -----
int NNODES = 20;        // Number of nodes per elements
int NELEM = 4;          // Number of elements

// ----- Helper Function Declarations -----

/**
 * @function C_P2
 * @brief Extraction operator matrix, P = 2
 * @param e : int : Element number
 * @returns C : Matrix3f : Extraction matrix
 */
Eigen::Matrix3f C_P2(int e) {
    Eigen::Matrix3f C;
    
    if (e == 1) {
        C(0,0) = 1.0;   C(0,1) = 0.0;   C(0,2) = 0.0;
        C(1,0) = 0.0;   C(1,1) = 1.0;   C(1,2) = 0.5;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 0.5;
    }
    else if (e == 4) {
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
    
    if (e == 1) {
        C(0,0) = 1.0;   C(0,1) = 0.0;   C(0,2) = 0.0;   C(0,3) = 0.0;
        C(1,0) = 0.0;   C(1,1) = 1.0;   C(1,2) = 0.5;   C(1,3) = 0.25;
        C(2,0) = 0.0;   C(2,1) = 0.0;   C(2,2) = 0.5;   C(2,3) = 7.0/12.0;
        C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = 1.0/6.0;
    }
    else if (e == 2) {
        C(0,0) = 0.25;     C(0,1) = 0.0;     C(0,2) = 0.0;     C(0,3) = 0.0;
        C(1,0) = 7.0/12.0; C(1,1) = 2.0/3.0; C(1,2) = 1.0/3.0; C(1,3) = 1.0/6.0;
        C(2,0) = 1.0/6.0;  C(2,1) = 1.0/3.0; C(2,2) = 2.0/3.0; C(2,3) = 2.0/3.0;
        C(3,0) = 0.0;      C(3,1) = 0.0;     C(3,2) = 0.0;     C(3,3) = 1.0/6.0;
    }
    else if (e == NELEM-1) {
        C(0,0) = 1.0/6.0; C(0,1) = 0.0;     C(0,2) = 0.0;     C(0,3) = 0.0;
        C(1,0) = 2.0/3.0; C(1,1) = 2.0/3.0; C(1,2) = 1.0/3.0; C(1,3) = 1.0/6.0;
        C(2,0) = 1.0/6.0; C(2,1) = 1.0/3.0; C(2,2) = 2.0/3.0; C(2,3) = 7.0/12.0;
        C(3,0) = 0.0;     C(3,1) = 0.0;     C(3,2) = 0.0;     C(3,3) = 1.0/4.0;
    }
    else if (e == NELEM) {
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

// ----- MAIN FUNCTION -----
/**
 * @function main
 * @brief Main driver for CE507_Coding2
 * @discussion Implements a BSpline basis function with extraction operators to Bezier curves.
 * @returns 1 : int : Program ran successfully
 * @returns 9 : int : Program terminated prematurely (user through debugging/errors)
 */
int main() {
    
    // BSpline Verification:
    //   The vectors of BSpline objects contain 3 and 4 instances of the BSpline basis
    //   functions, each with a different order. Their values are output into CSV files
    //   and verified before moving on.
    vector<BSpline> B_P2;
    B_P2.push_back(BSpline(2,1,NNODES));
    B_P2.push_back(BSpline(2,2,NNODES));
    B_P2.push_back(BSpline(2,3,NNODES));
    
    vector<BSpline> B_P3;
    B_P3.push_back(BSpline(3,1,NNODES));
    B_P3.push_back(BSpline(3,2,NNODES));
    B_P3.push_back(BSpline(3,3,NNODES));
    B_P3.push_back(BSpline(3,4,NNODES));
    
    ofstream domainOutput("domain.csv");
    ofstream B_P2Output1("B_P2_a1.csv");
    ofstream B_P2Output2("B_P2_a2.csv");
    ofstream B_P2Output3("B_P2_a3.csv");
    ofstream B_P3Output1("B_P3_a1.csv");
    ofstream B_P3Output2("B_P3_a2.csv");
    ofstream B_P3Output3("B_P3_a3.csv");
    ofstream B_P3Output4("B_P3_a4.csv");
    for (int i = 0; i < NNODES; i++) {
        domainOutput << B_P2[0].getDomain().getNodes()[i] << ',';
        B_P2Output1 << B_P2[0].getPoints()[i] << ',';
        B_P2Output2 << B_P2[1].getPoints()[i] << ',';
        B_P2Output3 << B_P2[2].getPoints()[i] << ',';
        B_P3Output1 << B_P3[0].getPoints()[i] << ',';
        B_P3Output2 << B_P3[1].getPoints()[i] << ',';
        B_P3Output3 << B_P3[2].getPoints()[i] << ',';
        B_P3Output4 << B_P3[3].getPoints()[i] << ',';
    }
    
    return 0;
}
