//
// Created by zz on 22/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_MATRIXUTILITY_H
#define IMPLEMENTATION_OPTIMIZATION_MATRIXUTILITY_H

#include <opencv2/opencv.hpp>
#include <vector>
using namespace cv;

class MatrixUtility {
public:
    static void matrixMultiply(const std::vector<double> &A, const unsigned int rowsA, const unsigned int colsA, const std::vector<double> &B, const unsigned int colsB, std::vector<double> &C) {
        Mat AMat(rowsA, colsA, CV_64F, const_cast<double*>(A.data()));
        Mat BMat(colsA, colsB, CV_64F, const_cast<double*>(B.data()));
        Mat CMat = AMat * BMat;
        for (int i = 0; i < rowsA; ++i) {
            for (int j = 0; j < colsB; ++j) {
                C[i*rowsA+j] = CMat.at<double>(i, j);
            }
        }
    }
    static double maximumNorm(const std::vector<double> vec) {
        double maxAbs = 0;
        for(auto val:vec) {
            if (fabs(val) > maxAbs) { maxAbs = fabs(val); }
        }
        return maxAbs;
    }
    static void applyOpenCVLeastSquaresSolver(const std::vector<double> &A, const std::vector<double> &negative_b, std::vector<double> &x) {
        cv::Mat bMat(negative_b.size(), 1, CV_64F, const_cast<double*>(negative_b.data()));
        cv::Mat AMat(negative_b.size(), x.size(), CV_64F, const_cast<double*>(A.data()));
        cv::Mat xMat(x.size(), 1, CV_64F);
        Mat AT = AMat.t();
        Mat ATA = AT*AMat;
        Mat ATb = AT*bMat;
        cv::solve(AMat, -bMat, xMat, cv::DECOMP_SVD);
        for (int i = 0; i < x.size(); ++i) {
            x[i] = xMat.at<double>(i, 0);
        }
    }

    static void test2() {
        Mat A1 = Mat::eye(10, 10, CV_32F);
        Mat C = A1.t(); // compute (A + lambda*I)^t * (A + lamda*I)
        std::cout<<"test2 done";
    }
    template <typename T>
    static void printVector(std::vector<T>& vec){
        std::cout<<"[";
        int i=0;
        for(auto ele: vec) {
            std::cout<<ele<<", ";
            i++;
        }
        std::cout<<"]";
    }
    static void printTriple(double* triple){
        std::cout<<"[";
        std::cout<<triple[0]<<", ";
        std::cout<<triple[1]<<", ";
        std::cout<<triple[2]<<", ";
        std::cout<<"]";
    }
};


#endif //IMPLEMENTATION_OPTIMIZATION_MATRIXUTILITY_H
