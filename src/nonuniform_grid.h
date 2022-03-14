//
// Created by Vlad Argunov on 13/11/2021.
//

#ifndef QUANTKIT_NONUNIFORM_GRID_H
#define QUANTKIT_NONUNIFORM_GRID_H

#include <Eigen/Dense>
#include <string>


void create_nonuniform_grid(Eigen::VectorXd & x_range, double x_min, double x_max, double alpha, double beta, int m);




#endif //QUANTKIT_NONUNIFORM_GRID_H
