//
// Created by zz on 21/04/2017.
//

#include "ExactLineSearch.h"
#include <cmath>

double ExactLineSearch::search() {
    if (derivativePhiAt0 >= 0) { return 0; }

    unsigned int k = 0;
    auto gamma = beta*derivativePhiAt0;
    double a = 0;
    double b = std::min(1.0, alpha_max);
    while ((phi(b)<=lambda(b)) && (derivativePhi(b)<=gamma)
           && (b<alpha_max && k<k_max)) {
        k += 1;
        a = b;
        b = std::min(2*b, alpha_max);
    }
    double alpha = b;
    double phiAlpha = phi(alpha);
    double derivativePhiAlpha = derivativePhi(alpha);
    double lambdaAlpha = lambda(alpha);
    while ((std::fabs(derivativePhiAlpha)>tau*std::fabs(derivativePhiAt0))
           && (b-a>epsilon) && (k<k_max)) {
        k += 1;
        refine(alpha, a, b);
        phiAlpha = phi(alpha);
        derivativePhiAlpha = derivativePhi(alpha);
        lambdaAlpha = lambda(alpha);
    }
    if (phi(alpha)>=phiAt0) { alpha = 0; }
    return alpha;
}

ExactLineSearch::ExactLineSearch(ResidualFunction &theFunction, double* x,
                                 std::vector<double>& h,
                                 const double alpha_max,
                                 const double rho,
                                 const double beta,
                               const unsigned short k_max,
                                 const double epsilon,
                                 const double tau):
        LineSearch(theFunction, x, h, alpha_max, rho, beta, k_max),
        epsilon(epsilon),
        tau(tau){}
