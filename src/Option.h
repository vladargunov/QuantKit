//
// Created by Vlad Argunov on 24/10/2021.
//

#ifndef QUANTKIT_OPTION_H
#define QUANTKIT_OPTION_H

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <algorithm>
#include <iomanip>


// Possible types of an Option
// Call vanilla - Call
// Put vanilla - Put


class Option {
private:
    int numdiff_t; // Number of steps in time
    int numdiff_s; // Number of steps in stock price

    bool fixed_step_time; // Specifies if the time step size is fixed
    bool fixed_step_stock; // Specifies if the stock step size is fixed

    double fixed_dt; // Change in time if the step is fixed
    double fixed_ds; // Change in price of stock if the step is fixed

    double s_max; // Maximum value_matrix of the stock
    double tn; // Time of expiration of an option

    double strike; // Strike of the option

    std::string type_option; // Type of the option considered
    bool european; // Specifies if the option is european = 1, or american = 0.

    bool stochastic_vol; // Specifies the type of volatility
    bool stochastic_rate; // Specifies the type of interest rate

    double deterministic_vol; // Value of the volatility if it is deterministic
    double deterministic_rate; // Value of the interest rate if it is deterministic

    std::string stock_boundary_condition; // Type of the boundary condition for a stock

    Eigen::MatrixXd solution_grid; // Solution matrix

    // Functions
    double parameter_a(int n);
    double parameter_b(int n);
    double parameter_c(int n);
    Eigen::MatrixXd implicit_matrix_vanilla();

    void create_boundary_condition_time(Eigen::VectorXd& boundary_vector);

    Eigen::MatrixXd solution_grid_vanilla(std::string method);

public:
    // Standard constructor
    explicit Option( std::string type_option_, bool european_ = true);
    // Set the parameters of the option
    void set_tn(double tn_);
    void set_s_max(double s_max_);
    void set_numdiff_t(int numdiff_t_);
    void set_numdiff_s(int numdiff_s_);
    void set_volatility(double vol_, bool stochastic_vol_ = 0 );
    void set_interest_rate(double rate_, bool stochastic_rate_ = 0);
    void set_strike(double strike_);
    void set_stock_boundary_condition(std::string stock_boundary_condition_);
    void fixed_difference_step(bool fixed_step_stock_ = 1, bool fixed_step_time_ = 1);

    void compute_solution_grid(std::string method="Implicit");
    void print_solution_grid(int precision = 3);


};


#endif //QUANTKIT_OPTION_H
