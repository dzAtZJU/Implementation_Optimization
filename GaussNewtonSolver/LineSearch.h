//
// Created by zz on 20/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_LINESEARCH_H
#define IMPLEMENTATION_OPTIMIZATION_LINESEARCH_H

#include "hw3_gn.h"
#include <vector>

class LineSearch {
public:
    LineSearch(ResidualFunction &theFunction, double* x, std::vector<double>& h, const double alpha_max=100, const double rho=0.001, const double beta=0.002,
                   const unsigned short k_max=1000);
    virtual double search() = 0;
protected:
    // Descent Function
    ResidualFunction& theFunction;
    double* x; /// constraint: array of size theFunction.nX()
    std::vector<double>& h; /// constraint: array of size theFunction.nX()
    // Parameter
    const double rho;/// constraint: (0,0.5)
    const double beta;/// constraint: (rho,1)
    const double alpha_max;/// constraint: [0,∞)
    const unsigned short k_max;
    // Cache
    std::vector<double> R;
    std::vector<double> J;
    double phiAt0;
    double derivativePhiAt0;
    double preAlpha;
    // Line Search Function
    double lambda(const double alpha) const;
    void eval(const double alpha);
    double phi(const double alpha);/// Expensive, use auxiliary variables to store
    double derivativePhi(const double alpha);/// Expensive, use auxiliary variables to store
    // Alg
    void refine(double &alpha, double &a, double &b);
};


#endif //IMPLEMENTATION_OPTIMIZATION_LINESEARCH_H
