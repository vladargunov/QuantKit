//
// Created by Vlad Argunov on 19/10/2021.
//

#include "Regression.h"

#include <Eigen/Dense>
#include <string>
#include <iostream>


Regression::Regression(std::string type) {
    type_regression = type;
}


void Regression::set_data(Eigen::MatrixXd x_values_, Eigen::VectorXd y_values_) {
    x_values = x_values_;
    y_values = y_values_;
}

void Regression::add_constant() {
    Eigen::MatrixXd x_values_with_constant;
    x_values_with_constant.resize(x_values.rows(), x_values.cols() + 1);
    x_values_with_constant.rightCols(x_values.cols()) = x_values;

    x_values_with_constant.col(0) = Eigen::VectorXd::Ones(x_values.rows());
    x_values = x_values_with_constant;

}

void Regression::print_data() {
    std::cout << "Data Matrix:\n";
    std::cout << x_values << "\n";
}

void Regression::print_responses() {
    std::cout << "Responses vector:\n";
    std::cout << y_values << "\n";
}

void Regression::solve() {
    if (type_regression == "Linear"){
        Eigen::MatrixXd computation_matrix;
        computation_matrix.resize(x_values.rows(), x_values.rows());
        computation_matrix = (x_values.transpose() * x_values).inverse();
        parameters = computation_matrix * x_values.transpose() * y_values;
    }
}

void Regression::print_parameters() {
    std::cout << "Parameters vector:\n";
    std::cout << parameters << "\n";
}
