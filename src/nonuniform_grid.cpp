//
// Created by Vlad Argunov on 13/11/2021.
//

#include "nonuniform_grid.h"

#include <Eigen/Dense>
#include <string>

void create_nonuniform_grid(Eigen::VectorXd & x_range, double x_min, double x_max, double alpha, double beta, int m) {
    /*The function transforms the grid in the non-uniform manner concentrating around the beta.
     * Parameter alpha is responsible for the fraction of points lying in the neighbourhood of the beta.
     * Important! The length of x_range must be of m + 1 exactly.*/

    double xi;
    double dxi = (asinh((x_max - beta)/alpha) - asinh((x_min - beta)/alpha)) / m;

    for (int i = 0; i < m + 1; ++i) {
        xi = asinh(- beta / alpha) + i * dxi;
        x_range(i,0) = beta + alpha * sinh(xi);
    }
}
