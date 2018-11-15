//   Created by Damyn Chipman on 11/15/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   Verifications.h
//   PROJECT:   CE507_Coding2

#ifndef Verifications_h
#define Verifications_h

#include <vector>
#include <Eigen/Dense>
#include "BSpline.h"
#include "DBSpline.h"

using namespace std;

void RunVerifications() {
    
    // BSpline Verification ===========================================================
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
    
    // Integration Points Verification ========================================================
    //    Compare values of BSpline vector to given values from Dr. Borden's calculations.
    //    Form of vectors:
    //      B(order 2, integration point[i]) = [# for basisID 1, # for basisID 2, # for basisID 3]
    //      B = [a = 1, a = 2, a = 3]
    //      DB(order 2, integration point[i]) = [# for basisID 1, # for basisID 2, # for basisID 3]
    //      DB = [a = 1, a = 2, a = 3]
    vector<float> intPoints = {static_cast<float>(-sqrt(3.0/5.0)), 0.0, static_cast<float>(sqrt(3.0/5.0))};
    BSpline B_P2_int(2, intPoints[0]);
    cout << "Integration Points Verification - B_P2:" << endl;
    cout << B_P2_int.getPoints()[0] << " ";
    cout << B_P2_int.getPoints()[1] << " ";
    cout << B_P2_int.getPoints()[2] << " " << endl << endl;
    
    cout << "Integration Points Verification - DB_P2:" << endl;
    DBSpline DB_P2_int(2, intPoints[0]);
    cout << DB_P2_int.getPoints()[0] << " ";
    cout << DB_P2_int.getPoints()[1] << " ";
    cout << DB_P2_int.getPoints()[2] << " " << endl << endl;
    
    // Extraction Operator Verification ========================================================
    //    Multiplying the BSpline vectors by the extraction operator matrix
    //    to obtain shape functions N. It's necessary to use Eigen's vector functionality
    //    in order to perform the multiplications. Checking with Dr. Borden's calculations.
    Eigen::VectorXf B_P2_vec = B_P2_int.getVector();
    Eigen::VectorXf N_P2_vec = C_P2(1)*B_P2_vec;
    cout << "Integration Points Verification - N_P2:" << endl;
    cout << N_P2_vec[0] << " ";
    cout << N_P2_vec[1] << " ";
    cout << N_P2_vec[2] << " " << endl << endl;
    
    Eigen::VectorXf DB_P2_vec = DB_P2_int.getVector();
    Eigen::VectorXf DN_P2_vec = C_P2(1)*DB_P2_vec;
    cout << "Integration Points Verification - DN_P2:" << endl;
    cout << DN_P2_vec[0] << " ";
    cout << DN_P2_vec[1] << " ";
    cout << DN_P2_vec[2] << " " << endl << endl;
}


#endif /* Verifications_h */
