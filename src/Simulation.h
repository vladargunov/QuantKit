//
// Created by Vlad Argunov on 03/01/2022.
//

#ifndef QUANTKIT_SIMULATION_H
#define QUANTKIT_SIMULATION_H

#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include <tuple>

class Simulation {

public:
    Eigen::VectorXd generate_uniform_variable(int n);
    Eigen::VectorXd generate_normal_variable(int n, double mu, double sigma);
    Eigen::MatrixXd generate_gbm(int n_processes, double tn, double stock_price, double mu, double sigma, int n_steps);
    std::tuple<Eigen::MatrixXd , Eigen::MatrixXd> generate_heston_process(int n_processes, double tn, double stock_price, double mu,
                                            double variance, double kappa, double theta, double sigma, double rho, int n_steps);

    Simulation() {};


};


#endif //QUANTKIT_SIMULATION_H
