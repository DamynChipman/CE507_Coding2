//   Created by Damyn Chipman on 10/26/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding2



#include <iostream>
#include <string>
#include <vector>
#include <fstream>
//#include <Eigen>
#include "Domain1D.h"
#include "BSpline.h"
//#include "FE1D.h"

using namespace std;
//using namespace Eigen;

int main(int argc, const char * argv[]) {
    
    int NNODES = 20; // Number of nodes per elements
    
    //BSpline B_P2(2, 1, NNODES);
    vector<BSpline> B_P2;
    B_P2.push_back(BSpline(2,1,NNODES));
    B_P2.push_back(BSpline(2,2,NNODES));
    B_P2.push_back(BSpline(2,3,NNODES));
    
    ofstream domainOutput("domain.csv");
    ofstream B_P2Output1("B_P2_a1.csv");
    ofstream B_P2Output2("B_P2_a2.csv");
    ofstream B_P2Output3("B_P2_a3.csv");
    for (int i = 0; i < NNODES; i++) {
        domainOutput << B_P2[0].getDomain().getNodes()[i] << ',';
        B_P2Output1 << B_P2[0].getPoints()[i] << ',';
        B_P2Output2 << B_P2[1].getPoints()[i] << ',';
        B_P2Output3 << B_P2[2].getPoints()[i] << ',';
    }
    
    return 0;
}
