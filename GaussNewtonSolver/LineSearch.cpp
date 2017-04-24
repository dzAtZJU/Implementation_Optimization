//
// Created by zz on 20/04/2017.
//

#include "LineSearch.h"
#include <numeric>
#include "MatrixUtility.h"
using namespace std;

LineSearch::LineSearch(ResidualFunction &theFunction, double* x, std::vector<double>& h, const double alpha_max, const double rho, const double beta,
                       const unsigned short k_max)
        : theFunction(theFunction), x(x), h(h), rho(rho), beta(beta), alpha_max(alpha_max), k_max(k_max) {
    R = std::vector<double>(theFunction.nR());
    J = std::vector<double>(theFunction.nX()*theFunction.nR());
    preAlpha = -1;
    phiAt0 = phi(0);
    derivativePhiAt0 = derivativePhi(0);
}

double LineSearch::lambda(const double alpha) const {
    return phiAt0 + rho*derivativePhiAt0*alpha;
}

double LineSearch::phi(const double alpha) {
    eval(alpha);
    return std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
}

double LineSearch::derivativePhi(const double alpha) {
    eval(alpha);
    vector<double> derivativeResidualFunction(theFunction.nX());
    MatrixUtility::matrixMultiply(R, 1, theFunction.nR(), J, theFunction.nX(), derivativeResidualFunction);
    for(auto& ele:derivativeResidualFunction) {
        ele *=2;
    }
    return inner_product(derivativeResidualFunction.begin(), derivativeResidualFunction.end(), h.begin(), 0.0);
}

void LineSearch::refine(double &alpha, double &a, double &b) {
    auto D = b - a;
    auto c = (phi(b) - phi(a) - D*derivativePhi(a))/(D*D);
    if (c>0) {
        alpha = a - derivativePhi(a)/(2*c);
        alpha = min(max(alpha, a+0.1*D),b-0.1*D);
    }
    else {
        alpha = (a+b)/2;
    }
    if (phi(alpha)<lambda(alpha)) {
        a = alpha;
    }
    else {
        b = alpha;
    }
}

void LineSearch::eval(const double alpha) {
    if (alpha == preAlpha) { return; }
    vector<double> nextX(theFunction.nX());
    for (int i = 0; i < nextX.size(); ++i) {
        nextX[i] = x[i] + alpha*h[i];
    }
    theFunction.eval(R.data(), J.data(), nextX.data());
    preAlpha = alpha;
}
