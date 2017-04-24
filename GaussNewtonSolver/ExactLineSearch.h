//
// Created by zz on 21/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_EXACTLINESEARCH_H
#define IMPLEMENTATION_OPTIMIZATION_EXACTLINESEARCH_H


#include "LineSearch.h"

class ExactLineSearch : public LineSearch {
public:
    double search() override;
    ExactLineSearch(ResidualFunction &theFunction, double* x, std::vector<double>& h, const double alpha_max=100, const double rho=0.001, const double beta=0.002,
    const unsigned short k_max=1000,const double epsilon=1e-5, const double tau=0.001);
private:
    const double epsilon;
    const double tau;
};


#endif //IMPLEMENTATION_OPTIMIZATION_EXACTLINESEARCH_H
