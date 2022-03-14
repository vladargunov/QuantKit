//
// Created by Vlad Argunov on 19/10/2021.
//

#ifndef QUANTKIT_REGRESSION_H
#define QUANTKIT_REGRESSION_H

#include <Eigen/Dense>
#include <string>
#include <iostream>

// Possible types of Regression
// Linear - Linear Regression


class Regression {
private:
    std::string type_regression;
    Eigen::MatrixXd x_values;
    Eigen::VectorXd y_values;
    Eigen::VectorXd parameters;
    bool constant;

public:
    // Standard constructor
    explicit Regression( std::string type);
    void set_data(Eigen::MatrixXd x_values_, Eigen::VectorXd y_values_);
    void add_constant();
    void print_data();
    void print_responses();
    void solve();
    void print_parameters();




};


#endif //QUANTKIT_REGRESSION_H
