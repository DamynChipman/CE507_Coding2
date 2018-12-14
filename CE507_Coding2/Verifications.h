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

vector<float> controlPoints(int n, float delX) {
    vector<float> CP(3);
    CP[0] = (n - 1)*delX;
    CP[1] = (delX/2.0)*(2*n + 1);
    CP[2] = n*delX;
    return CP;
}

void RunVerifications() {
    
    // BSpline Verification ===========================================================
    //   The vectors of BSpline objects contain 3 and 4 instances of the BSpline basis
    //   functions, each with a different order. Their values are output into CSV files
    //   and verified before moving on.
    vector<BSpline> B_P1;
    B_P1.push_back(BSpline(1,1,NNODES));
    B_P1.push_back(BSpline(1,2,NNODES));
    
    vector<BSpline> B_P2;
    B_P2.push_back(BSpline(2,1,NNODES));
    B_P2.push_back(BSpline(2,2,NNODES));
    B_P2.push_back(BSpline(2,3,NNODES));
    
    vector<BSpline> B_P3;
    B_P3.push_back(BSpline(3,1,NNODES));
    B_P3.push_back(BSpline(3,2,NNODES));
    B_P3.push_back(BSpline(3,3,NNODES));
    B_P3.push_back(BSpline(3,4,NNODES));
    
    vector<DBSpline> DB_P1;
    DB_P1.push_back(DBSpline(1,1,NNODES));
    DB_P1.push_back(DBSpline(1,2,NNODES));
    
    vector<DBSpline> DB_P2;
    DB_P2.push_back(DBSpline(2,1,NNODES));
    DB_P2.push_back(DBSpline(2,2,NNODES));
    DB_P2.push_back(DBSpline(2,3,NNODES));
    
    vector<DBSpline> DB_P3;
    DB_P3.push_back(DBSpline(3,1,NNODES));
    DB_P3.push_back(DBSpline(3,2,NNODES));
    DB_P3.push_back(DBSpline(3,3,NNODES));
    DB_P3.push_back(DBSpline(3,4,NNODES));
    
    ofstream domainOutput("domain.csv");
    ofstream B_P1Output1("B_P1_a1.csv");
    ofstream B_P1Output2("B_P1_a2.csv");
    ofstream B_P2Output1("B_P2_a1.csv");
    ofstream B_P2Output2("B_P2_a2.csv");
    ofstream B_P2Output3("B_P2_a3.csv");
    ofstream B_P3Output1("B_P3_a1.csv");
    ofstream B_P3Output2("B_P3_a2.csv");
    ofstream B_P3Output3("B_P3_a3.csv");
    ofstream B_P3Output4("B_P3_a4.csv");
    ofstream DB_P1Output1("DB_P1_a1.csv");
    ofstream DB_P1Output2("DB_P1_a2.csv");
    ofstream DB_P2Output1("DB_P2_a1.csv");
    ofstream DB_P2Output2("DB_P2_a2.csv");
    ofstream DB_P2Output3("DB_P2_a3.csv");
    ofstream DB_P3Output1("DB_P3_a1.csv");
    ofstream DB_P3Output2("DB_P3_a2.csv");
    ofstream DB_P3Output3("DB_P3_a3.csv");
    ofstream DB_P3Output4("DB_P3_a4.csv");
    for (int i = 0; i < NNODES; i++) {
        domainOutput << B_P2[0].getDomain().getNodes()[i] << ',';
        B_P1Output1 << B_P1[0].getPoints()[i] << ',';
        B_P1Output2 << B_P1[1].getPoints()[i] << ',';
        B_P2Output1 << B_P2[0].getPoints()[i] << ',';
        B_P2Output2 << B_P2[1].getPoints()[i] << ',';
        B_P2Output3 << B_P2[2].getPoints()[i] << ',';
        B_P3Output1 << B_P3[0].getPoints()[i] << ',';
        B_P3Output2 << B_P3[1].getPoints()[i] << ',';
        B_P3Output3 << B_P3[2].getPoints()[i] << ',';
        B_P3Output4 << B_P3[3].getPoints()[i] << ',';
        DB_P1Output1 << DB_P1[0].getPoints()[i] << ',';
        DB_P1Output2 << DB_P1[1].getPoints()[i] << ',';
        DB_P2Output1 << DB_P2[0].getPoints()[i] << ',';
        DB_P2Output2 << DB_P2[1].getPoints()[i] << ',';
        DB_P2Output3 << DB_P2[2].getPoints()[i] << ',';
        DB_P3Output1 << DB_P3[0].getPoints()[i] << ',';
        DB_P3Output2 << DB_P3[1].getPoints()[i] << ',';
        DB_P3Output3 << DB_P3[2].getPoints()[i] << ',';
        DB_P3Output4 << DB_P3[3].getPoints()[i] << ',';
    }
    
    // Local to Global Verification ===========================================================
    //    Create a mapping from the local domain in xi to the global domain in x.
    int p = 2;
    int NE = 10;
    float leftBound = 0.0;                                // Left domain bound
    float rightBound = 1.0;                               // Right domain bound
    float del_e = (rightBound - leftBound)/(NE);          // Element spacing
    int NINT = 3;                                         // Number of integration points for quadrature
    Eigen::VectorXf intPoints(NINT);                      // Vector of integration points
    intPoints(0) = -sqrt(3.0/5.0);                        // Integration point 1
    intPoints(1) = 0.0;                                   // Integration point 2
    intPoints(2) = sqrt(3.0/5.0);                         // Integration point 3
    int NKNOTS = (p+1)+(p+1)+(NE)-1;     
    Eigen::VectorXf knotVector(NKNOTS);                   // Knot vector
    int NDEL_E = 1;                                       // Spacing iterator
    for (int i = 0; i < NKNOTS; i++) {
        if (i < (p+1)) { knotVector(i) = leftBound; }
        if (i > ((p)+NE)) { knotVector(i) = rightBound; }
        if (i > p && i < (p+1+NE)) { knotVector(i) = NDEL_E*del_e; NDEL_E++; }
    }
    ofstream globalDomainOut("globalDomain.csv");
    ofstream N_P2_a1Out("N_P2_a1.csv");
    ofstream N_P2_a2Out("N_P2_a2.csv");
    ofstream N_P2_a3Out("N_P2_a3.csv");
    
    float x_xi = 0.0;
    for (int e = 0; e < NE; e++) {
        for (int i = 0; i < NINT; i++) {
            
            // Compute Bpi
            BSpline Bpi(p, intPoints(i));
            Eigen::VectorXf Bpi_vec = Bpi.getVector();
            
            // Compute Np
            Eigen::VectorXf Np_vec = C_P(p, e)*Bpi_vec;
            
            for (int A = 0; A < p+1; A++) {
                // Calc x_A
                float x_A = 0.0;
                for (int j = A+1; j < A + p + 1; j++) {
                    x_A = x_A + (1.0/p)*knotVector(j+e);
                }
                // Calc x
                x_xi = x_xi + x_A*Np_vec(A);
            }
            
            //cout << "x_xi = " << x_xi << " ";
            globalDomainOut << x_xi << ",";
            N_P2_a1Out << Np_vec(0) << ",";
            N_P2_a1Out << Np_vec(1) << ",";
            N_P2_a1Out << Np_vec(2) << ",";
        }
    }

    // Integration Points Verification ========================================================
    //    Compare values of BSpline vector to given values from Dr. Borden's calculations.
    //    Form of vectors:
    //      B(order 2, integration point[i]) = [# for basisID 1, # for basisID 2, # for basisID 3]
    //      B = [a = 1, a = 2, a = 3]
    //      DB(order 2, integration point[i]) = [# for basisID 1, # for basisID 2, # for basisID 3]
    //      DB = [a = 1, a = 2, a = 3]
    vector<float> intPoints2 = {static_cast<float>(-sqrt(3.0/5.0)), 0.0, static_cast<float>(sqrt(3.0/5.0))};
    BSpline B_P1_int1(1, intPoints[0]);
    BSpline B_P1_int2(1, intPoints[1]);
    BSpline B_P1_int3(1, intPoints[2]);
    cout << "Integration Points Verification - B_P1:" << endl;
    cout << "  Point 1: " << intPoints2[0] << endl;
    cout << B_P1_int1.getPoints()[0] << ", ";
    cout << B_P1_int1.getPoints()[1] << ", " << endl << endl;
    cout << "  Point 2: " << intPoints2[1] << endl;
    cout << B_P1_int2.getPoints()[0] << ", ";
    cout << B_P1_int2.getPoints()[1] << ", " << endl << endl;
    cout << "  Point 3: " << intPoints2[2] << endl;
    cout << B_P1_int3.getPoints()[0] << ", ";
    cout << B_P1_int3.getPoints()[1] << ", " << endl << endl;
    
    BSpline B_P2_int1(2, intPoints[0]);
    BSpline B_P2_int2(2, intPoints[1]);
    BSpline B_P2_int3(2, intPoints[2]);
    cout << "Integration Points Verification - B_P2:" << endl;
    cout << "  Point 1: " << intPoints2[0] << endl;
    cout << B_P2_int1.getPoints()[0] << ", ";
    cout << B_P2_int1.getPoints()[1] << ", ";
    cout << B_P2_int1.getPoints()[2] << ", " << endl << endl;
    cout << "  Point 2: " << intPoints2[1] << endl;
    cout << B_P2_int2.getPoints()[0] << ", ";
    cout << B_P2_int2.getPoints()[1] << ", ";
    cout << B_P2_int2.getPoints()[2] << ", " << endl << endl;
    cout << "  Point 3: " << intPoints2[2] << endl;
    cout << B_P2_int3.getPoints()[0] << ", ";
    cout << B_P2_int3.getPoints()[1] << ", ";
    cout << B_P2_int3.getPoints()[2] << ", " << endl << endl;
    
    BSpline B_P3_int1(3, intPoints[0]);
    BSpline B_P3_int2(3, intPoints[1]);
    BSpline B_P3_int3(3, intPoints[2]);
    cout << "Integration Points Verification - B_P2:" << endl;
    cout << "  Point 1: " << intPoints2[0] << endl;
    cout << B_P3_int1.getPoints()[0] << ", ";
    cout << B_P3_int1.getPoints()[1] << ", ";
    cout << B_P3_int1.getPoints()[2] << ", ";
    cout << B_P3_int1.getPoints()[3] << ", " << endl << endl;
    cout << "  Point 2: " << intPoints2[1] << endl;
    cout << B_P3_int2.getPoints()[0] << ", ";
    cout << B_P3_int2.getPoints()[1] << ", ";
    cout << B_P3_int2.getPoints()[2] << ", ";
    cout << B_P3_int2.getPoints()[3] << ", " << endl << endl;
    cout << "  Point 3: " << intPoints2[2] << endl;
    cout << B_P3_int3.getPoints()[0] << ", ";
    cout << B_P3_int3.getPoints()[1] << ", ";
    cout << B_P3_int3.getPoints()[2] << ", ";
    cout << B_P3_int3.getPoints()[3] << ", " << endl << endl;
    
    ofstream IP1_domainOut("IP_domain1.csv");
    ofstream IP2_domainOut("IP_domain2.csv");
    ofstream IP3_domainOut("IP_domain3.csv");
    ofstream IP_P1Out("IP_P1.csv");
    ofstream IP_P2Out("IP_P2.csv");
    ofstream IP_P3Out("IP_P3.csv");
    
    IP1_domainOut << intPoints2[0] << "," << intPoints2[0] << ",";
    IP1_domainOut << intPoints2[1] << "," << intPoints2[1] << ",";
    IP1_domainOut << intPoints2[2] << "," << intPoints2[2] << ",";
    
    IP2_domainOut << intPoints2[0] << "," << intPoints2[0] << "," << intPoints2[0] << ",";
    IP2_domainOut << intPoints2[1] << "," << intPoints2[1] << "," << intPoints2[1] << ",";
    IP2_domainOut << intPoints2[2] << "," << intPoints2[2] << "," << intPoints2[2] << ",";
    
    IP3_domainOut << intPoints2[0] << "," << intPoints2[0] << "," << intPoints2[0] << "," << intPoints2[0] << ",";
    IP3_domainOut << intPoints2[1] << "," << intPoints2[1] << "," << intPoints2[1] << "," << intPoints2[1] << ",";
    IP3_domainOut << intPoints2[2] << "," << intPoints2[2] << "," << intPoints2[2] << "," << intPoints2[2] << ",";
    
    IP_P1Out << B_P1_int1.getPoints()[0] << ",";
    IP_P1Out << B_P1_int1.getPoints()[1] << ",";
    IP_P1Out << B_P1_int2.getPoints()[0] << ",";
    IP_P1Out << B_P1_int2.getPoints()[1] << ",";
    IP_P1Out << B_P1_int3.getPoints()[0] << ",";
    IP_P1Out << B_P1_int3.getPoints()[1] << ",";
    
    IP_P2Out << B_P2_int1.getPoints()[0] << ",";
    IP_P2Out << B_P2_int1.getPoints()[1] << ",";
    IP_P2Out << B_P2_int1.getPoints()[2] << ",";
    IP_P2Out << B_P2_int2.getPoints()[0] << ",";
    IP_P2Out << B_P2_int2.getPoints()[1] << ",";
    IP_P2Out << B_P2_int2.getPoints()[2] << ",";
    IP_P2Out << B_P2_int3.getPoints()[0] << ",";
    IP_P2Out << B_P2_int3.getPoints()[1] << ",";
    IP_P2Out << B_P2_int3.getPoints()[2] << ",";
    
    IP_P3Out << B_P3_int1.getPoints()[0] << ",";
    IP_P3Out << B_P3_int1.getPoints()[1] << ",";
    IP_P3Out << B_P3_int1.getPoints()[2] << ",";
    IP_P3Out << B_P3_int1.getPoints()[3] << ",";
    IP_P3Out << B_P3_int2.getPoints()[0] << ",";
    IP_P3Out << B_P3_int2.getPoints()[1] << ",";
    IP_P3Out << B_P3_int2.getPoints()[2] << ",";
    IP_P3Out << B_P3_int2.getPoints()[3] << ",";
    IP_P3Out << B_P3_int3.getPoints()[0] << ",";
    IP_P3Out << B_P3_int3.getPoints()[1] << ",";
    IP_P3Out << B_P3_int3.getPoints()[2] << ",";
    IP_P3Out << B_P3_int3.getPoints()[3] << ",";
    
    
    // Extraction Operator Verification ========================================================
    //    Multiplying the BSpline vectors by the extraction operator matrix
    //    to obtain shape functions N. It's necessary to use Eigen's vector functionality
    //    in order to perform the multiplications. Checking with Dr. Borden's calculations.
//    Eigen::VectorXf B_P2_vec = B_P2_int.getVector();
//    Eigen::VectorXf N_P2_vec = C_P2(1)*B_P2_vec;
//    cout << "Integration Points Verification - N_P2:" << endl;
//    cout << N_P2_vec[0] << " ";
//    cout << N_P2_vec[1] << " ";
//    cout << N_P2_vec[2] << " " << endl << endl;
//
//    Eigen::VectorXf DB_P2_vec = DB_P2_int.getVector();
//    Eigen::VectorXf DN_P2_vec = C_P2(1)*DB_P2_vec;
//    cout << "Integration Points Verification - DN_P2:" << endl;
//    cout << DN_P2_vec[0] << " ";
//    cout << DN_P2_vec[1] << " ";
//    cout << DN_P2_vec[2] << " " << endl << endl;
}


#endif /* Verifications_h */
