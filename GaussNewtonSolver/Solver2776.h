//
// Created by zz on 21/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_SOLVER2776_H
#define IMPLEMENTATION_OPTIMIZATION_SOLVER2776_H

#include "hw3_gn.h"
#include <vector>

//Implementation of GaussNewtonSolver
class Solver2776: GaussNewtonSolver {
public:
    double solve(
            ResidualFunction *f, // 目标函数
            double *X,           // 输入作为初值，输出作为结果
            GaussNewtonParams param = GaussNewtonParams(), // 优化参数
            GaussNewtonReport *report = nullptr // 优化结果报告
    ) override;
};


#endif //IMPLEMENTATION_OPTIMIZATION_SOLVER2776_H
