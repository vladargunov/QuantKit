//
// Created by Vlad Argunov on 17/01/2022.
//

#ifndef QUANTKIT_OPTIMALEXECUTION_H
#define QUANTKIT_OPTIMALEXECUTION_H

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "Simulation.h"

// Impact: Linear, Nonlinear
// Total time: 1
// Initial quantity: 1
class OptimalExecution {
private:
    int num_diff_t;
    int num_diff_q;

    Eigen::VectorXd stock;

    Eigen::VectorXd cash;

    Eigen::VectorXd inventory;

    // Value function at any time and quantity
    Eigen::MatrixXd value_matrix;
    // Value function along the executed directory
    Eigen::VectorXd value_full;

    // Matrix of optimal trading speed that captures the speed at any quantity
    Eigen::MatrixXd optimal_speed_matrix;
    // Vectors of optimal trading speed taken along the executed directory, i.e. according to the chosen inventory path
    Eigen::VectorXd optimal_speed;

    Simulation sim;

    double a_coeff;
    double b_coeff;
    double k_coeff;
    double sigma_coeff;
    double phi_coeff;

public:
    OptimalExecution(double a_coeff_, double b_coeff_, double sigma_coeff_, double k_coeff_, double phi_coeff_);
    void compute_liquidation(const std::string& trading_speed, double impact_nonlinearity = 1);
    void set_num_steps(int num_diff_t_, int num_diff_q_);

    Eigen::VectorXd get_value_process() {return value_full;};
    Eigen::MatrixXd get_value_matrix() {return value_matrix;};
    Eigen::VectorXd get_optimal_speed_process() {return optimal_speed;};
    Eigen::MatrixXd get_optimal_speed_matrix() {return optimal_speed_matrix;};
    Eigen::VectorXd get_stock_process() {return stock;};
    Eigen::VectorXd get_cash_process() {return cash;};
    Eigen::VectorXd get_inventory_process() {return inventory;};




};


#endif //QUANTKIT_OPTIMALEXECUTION_H
