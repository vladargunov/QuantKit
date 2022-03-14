//
// Created by Vlad Argunov on 03/01/2022.
//

#include "Simulation.h"

#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <tuple>

Eigen::VectorXd Simulation::generate_uniform_variable(int n) {
        std::random_device rd;
        std::mt19937 gen(rd());
        Eigen::VectorXd uniform_numbers;
        uniform_numbers.resize(n);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        for (int i = 0; i < n; ++i) {
            uniform_numbers(i) = uniform(gen);
        }
        return uniform_numbers;
}

Eigen::VectorXd Simulation::generate_normal_variable(int n, double mu, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    Eigen::VectorXd normal_numbers;
    normal_numbers.resize(n);
    std::normal_distribution<> normal{mu,sigma};
    for (int i = 0; i < n; ++i) {
        normal_numbers(i) = normal(gen);
    }
    return normal_numbers;
}

Eigen::MatrixXd Simulation::generate_gbm(int n_processes, double tn, double stock_price, double mu, double sigma, int n_steps) {
//    Generates the matrix of Geometric brownian motions for each row and a time step for each column
    double time_step = tn / n_steps;
    Eigen::MatrixXd simulated_processes;
    Eigen::Index n_rows = n_processes;
    Eigen::Index n_cols = n_steps + 1;
    simulated_processes.resize(n_rows, n_cols);
    simulated_processes.leftCols(1) = Eigen::MatrixXd::Constant(n_rows, 1, stock_price);

    Eigen::VectorXd normal_rv;
    normal_rv.resize(n_processes);
    for (int col = 0; col < n_cols - 1; ++col) {
        normal_rv = generate_normal_variable(n_processes, 0, 1);
        for (int row = 0; row < n_rows; ++row) {
            simulated_processes(row, col + 1) = simulated_processes(row, col)
                    * exp(mu - 0.5 * sigma * sigma * time_step + sigma * sqrt(time_step) * normal_rv(row) );
        }

    }

    return simulated_processes;
}

std::tuple<Eigen::MatrixXd , Eigen::MatrixXd>
Simulation::generate_heston_process(int n_processes, double tn, double stock_price, double mu, double variance,
                                    double kappa, double theta, double sigma, double rho, int n_steps) {
    //    Generates the matrix of Heston processes for each row and a time step for each column
    // Here variance is the starting value_matrix of the volatility, theta is the long-term level,
    // kappa is the speed of convergence, and sigma is the volatility of volatility
    double time_step = tn / n_steps;
    Eigen::MatrixXd simulated_stock_processes, simulated_var_process;
    Eigen::Index n_rows = n_processes;
    Eigen::Index n_cols = n_steps + 1;

    simulated_stock_processes.resize(n_rows, n_cols);
    simulated_var_process.resize(n_rows, n_cols);

    simulated_stock_processes.leftCols(1) = Eigen::MatrixXd::Constant(n_rows, 1, stock_price);
    simulated_var_process.leftCols(1) = Eigen::MatrixXd::Constant(n_rows, 1, variance);

    Eigen::VectorXd normal_rv1, normal_rv2;
    normal_rv1.resize(n_processes);
    normal_rv2.resize(n_processes);
    for (int col = 0; col < n_cols - 1; ++col) {
        normal_rv1 = generate_normal_variable(n_processes, 0, 1);
        for (int row = 0; row < n_rows; ++row) {
            simulated_var_process(row, col + 1) = simulated_var_process(row, col) +
                                                  kappa * (theta - simulated_var_process(row, col)) * time_step
                                                  + sigma * sqrt(simulated_var_process(row, col) * time_step)
                                                  * (rho * normal_rv1(row) + sqrt(1 - rho * rho) * normal_rv2(row));


            simulated_stock_processes(row, col + 1) = simulated_stock_processes(row, col) * (1 + mu * time_step)
                    + simulated_stock_processes(row, col) * sqrt(simulated_var_process(row, col) * time_step) * normal_rv1(row);
        }

    }
    return {simulated_stock_processes, simulated_var_process};
}
