//   Created by Damyn Chipman on 10/26/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   Domain1D.h
//   PROJECT:   CE507_Coding2

#ifndef Domain1D_h
#define Domain1D_h

#include <sstream>

class Domain1D {
    
public:
    
    Domain1D(float leftBound, float rightBound, int n, std::string gridType) : lBound_(leftBound), rBound_(rightBound), N_(n), type_(gridType) {
        
        nodes_ = new float[N_];
        
        // TODO: Use correct spacing
        if (type_ == "edge") {
            del_x_ = (rBound_ - lBound_)/(N_ - 1);
        }
        else if (type_ == "centered") {
            del_x_ = (rBound_ - lBound_)/(N_ - 1);
        }
        else if (type_ == "centered ghost") {
            del_x_ = (rBound_ - lBound_)/(N_ - 1);
        }
        
        for (int i = 0; i < N_; i++) {
            nodes_[i] = i*del_x_ + leftBound;
        }
    }
    
    Domain1D() {}
    
    //
    void chageGridType(std::string newGrid) {
        // TODO: Build this function
    }
    
    //
    std::string viewNodes() {
        std::stringstream toOutput;
        for (int i = 0; i < N_; i++) {
            toOutput << "Node i: " << i << " = " << nodes_[i] << std::endl;
        }
        return toOutput.str();
    }
    
    // ----- Accessor Functions -----
    float getLBound() const { return lBound_; }
    float getRBound() const { return rBound_; }
    int getN() const { return N_; }
    float* getNodes() const { return nodes_; }
    std::string getType() const { return type_; }
    float getDelX() const { return del_x_; }
    
private:
    float lBound_;
    float rBound_;
    int N_;
    float* nodes_;
    std::string type_;
    float del_x_;
};

#endif /* Domain_h */
