//
// Created by zz on 21/04/2017.
//

#include <vector>
#include "Solver2776.h"
#include <numeric>
#include "SoftLineSearch.h"
#include "ExactLineSearch.h"
#include "MatrixUtility.h"

double Solver2776::solve(ResidualFunction *f, double *X, GaussNewtonParams param, GaussNewtonReport *report) {
    auto nR = f->nR();
    auto nX = f->nX();
    // 分配好 R 和 J需要的空间
    std::vector<double> R(nR);
    std::vector<double> J(nR*nX);
    //J * DeltaX = -R
    std::vector<double> DeltaX(nX);
    //迭代次数
    unsigned int n = 0;
    while (n < param.max_iter) {
        double x = X[0], y = X[1], z = X[2];
        //计算 Residual Function的余项和雅克比
        f->eval(R.data(), J.data(), X);
        //如果余项误差很小了，完成
        if (MatrixUtility::maximumNorm(R)<=param.residual_tolerance) {
            if(report!= nullptr) {
                report->stop_type = GaussNewtonReport::STOP_RESIDUAL_TOL;
            }
            break;
        }
        MatrixUtility::applyOpenCVLeastSquaresSolver(J, R, DeltaX);
        //如果下降方向很平缓了，完成
        if (MatrixUtility::maximumNorm(DeltaX)<=param.gradient_tolerance) {
            if(report!= nullptr) {
                report->stop_type = GaussNewtonReport::STOP_GRAD_TOL;
            }
            break;
        }
        //根据输入的参数，选择SoftLineSearch或者ExactLineSearch，并调用search()计算步长alpha
        double alpha = param.exact_line_search? ExactLineSearch(*f, X, DeltaX).search():
                       SoftLineSearch(*f, X, DeltaX).search();
        //按照下降方向和步长下降后的新的X
        for (int i = 0; i < nX; ++i) {
            X[i] += alpha*DeltaX[i];
        }
        if(param.verbose) {
            std::cout<<n<<"  ";
            std::cout<<"alpha: "<<alpha<<"  ";
            std::cout<<"deltaX: "; MatrixUtility::printVector(DeltaX); std::cout<<"  ";
            std::cout<<"movedX: "; MatrixUtility::printTriple(X); std::cout<<"  ";
            std::cout<<std::endl;
        }
        n++;
    }

    if(report!= nullptr) {
        if (n >= param.max_iter) {
            report->stop_type = GaussNewtonReport::STOP_NO_CONVERGE;
        }
        report->n_iter = n;
    }
    //计算该函数的最小值并返回
    return std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
}