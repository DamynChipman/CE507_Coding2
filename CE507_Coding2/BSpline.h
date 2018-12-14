//   Created by Damyn Chipman on 10/31/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   BSpline.h
//   PROJECT:   CE507_Coding2

#ifndef BSpline_h
#define BSpline_h

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "Domain1D.h"

class BSpline {
    
public:
    
    BSpline(int p, int a, int N) : order_(p), basisID_(a), domain_(Domain1D(-1.0, 1.0, N, "edge")), N_(N) {
        
        // Set up points within domain of BSpline
        points_ = new float[domain_.getN()];
        for (int i = 0; i < domain_.getN(); i++) {
            points_[i] = eval(domain_.getNodes()[i]);
        }
    };
    
    BSpline(int p, int a, int NINT, std::vector<float> intPoints) : order_(p), basisID_(a), domain_(Domain1D(-1.0, 1.0, NINT, "edge")), N_(NINT) {
        
        points_ = new float[NINT];
        for (int i = 0; i < intPoints.size(); i++) {
            points_[i] = eval(intPoints[i]);
        }
    }
    
    BSpline(int p, float intPoint) : order_(p), N_(p+1) {
        
        points_ = new float[p+1];
        for (int i = 0; i < p+1; i++) {
            basisID_ = i+1;
            points_[i] = eval(intPoint);
        }
    }
    
    float eval(float X) {
        float res = (1.0/(2.0*order_))*(factorial(order_)/(factorial(basisID_ - 1.0) * factorial(order_ + 1.0 - basisID_))) * pow(1 - X, order_ - (basisID_ - 1)) * pow(1 + X, basisID_ - 1);
        
        if (order_ == 3) {
            if (basisID_ == 1) {
                res = (3.0/4.0)*(1.0/6.0)*pow(1 - X,3);
            }
            else if (basisID_ == 2) {
                res = (3.0/4.0)*(1.0/2.0)*(1 + X)*pow(1 - X,2);
            }
            else if (basisID_ == 3) {
                res = (3.0/4.0)*(1.0/2.0)*(1 - X)*pow(1 + X,2);
            }
            else if (basisID_ == 4) {
                res = (3.0/4.0)*(1.0/6.0)*pow(1 + X,3);
            }
            else {
                
            }
        }
        return res;
    }
    
    float* operator*(float a) {
        float* newPoints = points_;
        for (int i = 0; i < domain_.getN(); i++) {
            newPoints[i] = points_[i]*a;
        }
        return newPoints;
    }
    
    // Accessor Functions
    float* getPoints() { return points_; }
    Domain1D getDomain() { return domain_; }
    int getBasisID() { return basisID_; }
    int getOrder() { return order_; }
    Eigen::VectorXf getVector() {
        Eigen::VectorXf vecPoints(N_);
        for (int i = 0; i < N_; i++) {
            vecPoints(i) = points_[i];
        }
        return vecPoints;
    }
    
    ~BSpline() {};
    
private:
    
    float factorial(float n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1)*n; }
    
    float* points_;
    int basisID_ = 0;
    int order_;
    int N_;
    Domain1D domain_;
    
};


#endif /* BSpline_h */
