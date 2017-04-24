//
// Created by zz on 24/04/2017.
//

#include "GaussNewtonSolver/EllipsoidFitting.h"
#include "GaussNewtonSolver/Solver2776.h"
using namespace std;

//Gauss Newton Solver Use Example
int main() {
    EllipsoidFitting ellipsoidFitting("/Users/tgbus/ClionProjects/Implementation_Optimization/GaussNewtonSolver/ellipse753.txt");
    Solver2776 solver;
    vector<double> X{1,1,1};
    auto gnParams = GaussNewtonParams(); gnParams.exact_line_search = true; gnParams.verbose = true;
    auto report = GaussNewtonReport();
    auto minimum = solver.solve(&ellipsoidFitting, X.data(), gnParams, &report);
    return 0;
}