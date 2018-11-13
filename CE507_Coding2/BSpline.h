//   Created by Damyn Chipman on 10/31/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   BSpline.h
//   PROJECT:   CE507_Coding2

#ifndef BSpline_h
#define BSpline_h

#include <cmath>
#include "Domain1D.h"

class BSpline {
    
public:
    
    BSpline(int p, int a, int N) : order_(p), basisID_(a), domain_(Domain1D(-1.0, 1.0, N, "edge")) {
        
        // Set up points within domain of BSpline
        points_ = new float[domain_.getN()];
        for (int i = 0; i < domain_.getN(); i++) {
            points_[i] = eval(domain_.getNodes()[i]);
        }
    };
    
    float eval(float X) {
        float res = (1.0/(2.0*order_))*(factorial(order_)/(factorial(basisID_ - 1.0) * factorial(order_ + 1.0 - basisID_)));
        res = res * pow(1 - X, order_ - (basisID_ - 1)) * pow(1 + X, basisID_ - 1);
        return res;
    }
    
    
    // Accessor Functions
    float* getPoints() { return points_; }
    Domain1D getDomain() { return domain_; }
    
    ~BSpline() {};
    
private:
    
    float factorial(float n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1)*n; }
    
    float* points_;
    int basisID_ = 0;
    int order_;
    Domain1D domain_;
    
};


#endif /* BSpline_h */
