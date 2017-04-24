//
// Created by zz on 21/04/2017.
//

#ifndef IMPLEMENTATION_OPTIMIZATION_EllipsoidFitting_H
#define IMPLEMENTATION_OPTIMIZATION_EllipsoidFitting_H

#include "hw3_gn.h"
#include <vector>
#include <string>

class EllipsoidFitting: public ResidualFunction {
public:
    EllipsoidFitting(std::string datasFile);
    int nX() const override;
    int nR() const override;
    void eval(double *R, double *J, double *X) override;
    void readInData(std::string datasFile);
private:
    std::vector<std::vector<double>> xyz_triples;
    double ellipsoid(std::vector<double> XYZ, const double* const ABC) const;
    double aSquareDivideBybCubic(double a, double b) const;
    double test1() const;
};


#endif //IMPLEMENTATION_OPTIMIZATION_EllipsoidFitting_H
