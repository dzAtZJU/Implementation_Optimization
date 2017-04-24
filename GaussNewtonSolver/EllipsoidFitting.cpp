//
// Created by zz on 21/04/2017.
//

#include "EllipsoidFitting.h"
#include <fstream>
#include <assert.h>

EllipsoidFitting::EllipsoidFitting(std::string datasFile) {
    readInData(datasFile);
}

void EllipsoidFitting::readInData(std::string datasFile) {
    std::ifstream dataFile;
    dataFile.open(datasFile, std::ios::in);
    if(!dataFile.is_open()) { assert(true); }
    double x, y, z;
    while(dataFile>>x, dataFile>>y, dataFile>>z) {
        xyz_triples.push_back(std::vector<double>{x, y, z});
        assert(xyz_triples.back().size()==3);
    }
    dataFile.close();
}

int EllipsoidFitting::nX() const {
    return 3;
}

int EllipsoidFitting::nR() const {
    return xyz_triples.size();
}

void EllipsoidFitting::eval(double *R, double *J, double *X) {
    for (int i = 0; i < nR(); ++i) {
        auto xyz_triple = xyz_triples[i];
        R[i] = ellipsoid(xyz_triples[i], X);
        for (int j = 0; j < nX(); ++j) {
            auto axis = xyz_triples[i][j];
            auto para = X[j];
            auto der = -2*aSquareDivideBybCubic(axis, para);
            auto index = i*nX()+j;
            J[index] = der;
        }
        test1();
    }
}

double EllipsoidFitting::ellipsoid(std::vector<double> XYZ, const double* const ABC) const {
    auto A = ABC[0];
    auto B = ABC[1];
    auto C = ABC[2];
    auto x = XYZ[0];
    auto y = XYZ[1];
    auto z = XYZ[2];
    return x*x/(A*A) + y*y/(B*B) + z*z/(C*C) - 1;
}

double EllipsoidFitting::aSquareDivideBybCubic(double a, double b) const {
    return a*a/(b*b*b);
}

double EllipsoidFitting::test1() const {
    for(auto& ele: xyz_triples) {
        assert(ele.size()==3);
    }
    return 0 ;
}

