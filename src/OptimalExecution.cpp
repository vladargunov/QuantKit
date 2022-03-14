//
// Created by Vlad Argunov on 17/01/2022.
//

#include "OptimalExecution.h"

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "Simulation.h"

OptimalExecution::OptimalExecution(double a_coeff_, double b_coeff_, double sigma_coeff_, double k_coeff_, double phi_coeff_) {
    a_coeff = a_coeff_;
    b_coeff = b_coeff_;
    sigma_coeff = sigma_coeff_;
    k_coeff = k_coeff_;
    phi_coeff = phi_coeff_;
}


void OptimalExecution::compute_liquidation(const std::string& trading_speed, double impact_nonlinearity) {
    int num_rows = num_diff_q + 1;
    int num_cols = num_diff_t + 1;
    sim = Simulation();
    value_matrix.resize(num_rows, num_cols);
    optimal_speed_matrix.resize(num_rows, num_cols);
    double diff_q = 1.0 / num_diff_q;
    double diff_t = 1.0 / num_diff_t;
    double quantity;
    double h_diff_q_up;
    double h_diff_q_down;
    double partial_diff_q;

    Eigen::VectorXd diffusion; diffusion.resize(num_cols - 1);
    diffusion = sim.generate_normal_variable(num_cols - 1, 0, pow(diff_t, 2));

    stock.resize(num_cols); inventory.resize(num_cols);
    cash.resize(num_cols); optimal_speed.resize(num_cols);
    value_full.resize(num_cols);

    if (trading_speed == "Nonlinear" || trading_speed == "Linear") {

        // That part computes the matrix of value_matrix function for each value of time and quantity.
        for (int col = (num_cols - 1); col > -1; --col) { // Time space
            for (int row = 0; row < num_rows; ++row) { // Stock space
                quantity = row * diff_q;

                if (trading_speed == "Linear"){
                    double gamma = sqrt(phi_coeff / k_coeff);
                    double zeta = (a_coeff - 0.5 * b_coeff + sqrt(k_coeff * phi_coeff)) / (a_coeff - 0.5 * b_coeff - sqrt(k_coeff * phi_coeff));
                    double remaining_time = (num_cols - col) * diff_t;
                    optimal_speed_matrix(row, col) = gamma * (zeta * exp(gamma * remaining_time) + exp(- gamma * remaining_time)) / (zeta * exp(gamma) - exp(- gamma));
                }

                if (col == (num_cols - 1)){
                    value_matrix(row, col) = - a_coeff * quantity * quantity;

                    if (trading_speed == "Nonlinear") {
                        optimal_speed_matrix(row, col) = pow(abs(-(b_coeff * quantity - 2 * a_coeff * quantity) / ((1 + a_coeff) * k_coeff))
                                ,1.0 / impact_nonlinearity);
                    }
                } else {

                    if (row == 0){
                        h_diff_q_up =  value_matrix(row + 1, col + 1);
                        h_diff_q_down =  value_matrix(row, col + 1);
                        partial_diff_q = (h_diff_q_up - h_diff_q_down) / (diff_q);
                    } else if (row == (num_rows - 1)) {
                        h_diff_q_up =  value_matrix(row, col + 1);
                        h_diff_q_down =  value_matrix(row - 1, col + 1);
                        partial_diff_q = (h_diff_q_up - h_diff_q_down) / (diff_q);
                    } else {
                        h_diff_q_up =  value_matrix(row + 1, col + 1);
                        h_diff_q_down =  value_matrix(row - 1, col + 1);
                        partial_diff_q = (h_diff_q_up - h_diff_q_down) / (2 * diff_q);
                    }

                    value_matrix(row, col) = value_matrix(row, col + 1) -
                            diff_t * (phi_coeff * pow(quantity, 2) - a_coeff * k_coeff *
                            pow(abs(-(b_coeff * quantity + partial_diff_q ) / ((1 + a_coeff) * k_coeff))
                                ,1.0 + 1.0/ impact_nonlinearity));


                    if (trading_speed == "Nonlinear"){
                        optimal_speed_matrix(row, col) = pow(abs(-(b_coeff * quantity + partial_diff_q) / ((1 + a_coeff) * k_coeff))
                                ,1.0 / impact_nonlinearity);
                    }

                }
            }
        }

        // That part simulates stock, cash, and inventory processes
        for (int step = 0; step < num_cols; ++step) {
            if (step == 0){
                inventory(step) = 1;
                stock(step) = 1;
                cash(step) = 0;

                optimal_speed(step) = optimal_speed_matrix(num_rows - 1, step);
                value_full(step) = cash(step) + inventory(step) * stock(step) + value_matrix(num_rows - 1, step);


            } else {
                inventory(step) = std::max(inventory(step - 1) - optimal_speed(step - 1) * diff_t, 0.0);
                stock(step) = stock(step - 1) - b_coeff * optimal_speed(step - 1) * diff_t
                        + sigma_coeff * diffusion(step - 1);
                cash(step) = cash(step - 1) + (stock(step - 1)
                        - k_coeff * pow(optimal_speed(step - 1), impact_nonlinearity)) * optimal_speed(step - 1) * diff_t;

                // Determine the relevant trading speed and value function
                for (int row = 0; row < num_rows; ++row) {
                    if ((row * diff_q <= inventory(step)) && ((row + 1) * diff_q > inventory(step))){
                        if (row != (num_rows - 1)){
                            optimal_speed(step) =
                                    (((row + 1) * diff_q - inventory(step)) / diff_q) * optimal_speed_matrix(row, step) +
                                    ((inventory(step) - row * diff_q) / diff_q) * optimal_speed_matrix(row + 1, step);


                            value_full(step) = cash(step) + inventory(step) * stock(step) +
                                    (((row + 1) * diff_q - inventory(step)) / diff_q) * value_matrix(row, step) +
                                    ((inventory(step) - row * diff_q) / diff_q) * value_matrix(row + 1, step);
                        } else {
                            optimal_speed(step) = optimal_speed_matrix(num_rows - 1, step);
                            value_full(step) = cash(step) + inventory(step) * stock(step) + value_matrix(num_rows - 1, step);
                        }

                    }
                }
            }
        }
    }

    if (inventory(num_cols - 1) > 0.0){
        cash(num_cols - 1) += inventory(num_cols - 1) * (stock(num_cols - 1) - a_coeff * inventory(num_cols - 1));
    }
    }





void OptimalExecution::set_num_steps(int num_diff_t_, int num_diff_q_) {
    num_diff_t = num_diff_t_;
    num_diff_q = num_diff_q_;
}
