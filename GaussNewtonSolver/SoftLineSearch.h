//
// Created by zz on 21/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_SOFTLINESEARCH_H
#define IMPLEMENTATION_OPTIMIZATION_SOFTLINESEARCH_H


#include "LineSearch.h"

class SoftLineSearch: public LineSearch {
public:
    SoftLineSearch(ResidualFunction &theFunction, double* x, std::vector<double>& h, const double alpha_max=100, const double rho=0.001, const double beta=0.002,
    const unsigned short k_max=1000);
    double search()  override;
};

#endif //IMPLEMENTATION_OPTIMIZATION_SOFTLINESEARCH_H
