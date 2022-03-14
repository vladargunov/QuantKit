//
// Created by Vlad Argunov on 24/10/2021.
//

#include "Option.h"
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <algorithm>
#include <iomanip>


Option::Option(std::string type_option_, bool european_){
    type_option = type_option_;
    european = european_;
}

void Option::set_tn(double tn_) {
    tn = tn_;

}

void Option::set_s_max(double s_max_) {
    s_max = s_max_;
}

void Option::set_numdiff_t(int numdiff_t_) {
    numdiff_t = numdiff_t_;
}

void Option::set_numdiff_s(int numdiff_s_) {
    numdiff_s = numdiff_s_;
    }

void Option::set_volatility(double vol_, bool stochastic_vol_) {
    stochastic_vol = stochastic_vol_;
    if (stochastic_vol == 0){
        deterministic_vol = vol_;
    }
    // Code for implementation of stochastic volatility
}

void Option::set_interest_rate(double rate_, bool stochastic_rate_) {
    stochastic_rate = stochastic_rate_;
    if (stochastic_rate == 0){
        deterministic_rate = rate_;
    }
    // Code for implementation of stochastic interest rate
}

void Option::set_strike(double strike_) {
    strike = strike_;
}

void Option::set_stock_boundary_condition(std::string stock_boundary_condition_) {
    stock_boundary_condition = stock_boundary_condition_;
}


double Option::parameter_a(int n) {
    if (stochastic_vol == 0 && stochastic_rate == 0 && fixed_step_time == 1){
        return 0.5 * (deterministic_rate * n - deterministic_vol * deterministic_vol * n * n) * fixed_dt;
    }

}

double Option::parameter_b(int n) {
    if (stochastic_vol == 0 && stochastic_rate == 0 && fixed_step_time == 1) {
        return (deterministic_vol * deterministic_vol * n * n + deterministic_rate) * fixed_dt;
    }
}

double Option::parameter_c(int n) {
    if (stochastic_vol == 0 && stochastic_rate == 0 && fixed_step_time == 1) {
        return -0.5 * (deterministic_vol * deterministic_vol * n * n + deterministic_rate * n) * fixed_dt;
    }
}


Eigen::MatrixXd Option::implicit_matrix_vanilla() {
    double a;
    double b;
    double c;
    Eigen::MatrixXd computation_matrix;
    computation_matrix.resize(numdiff_s - 1, numdiff_s - 1);
    for (int row = 0; row < computation_matrix.rows(); ++row) {
        if (row == 0){
            b = parameter_b(row + 1);
            c = parameter_c(row + 1);
            computation_matrix(row,row) = 1 + b;
            computation_matrix(row,row + 1) = c;
        } else if (row == computation_matrix.rows() - 1) {
            b = parameter_b(row + 1);
            a = parameter_a(row + 1);
            computation_matrix(row, row) = 1 + b;
            computation_matrix(row, row - 1) = a;
        } else {
            a = parameter_a(row + 1);
            b = parameter_b(row + 1);
            c = parameter_c(row + 1);
            computation_matrix(row, row) = 1 + b;
            computation_matrix(row, row - 1) = a;
            computation_matrix(row, row + 1) = c;
        }
    }
    return computation_matrix;
}

void Option::fixed_difference_step(bool fixed_step_stock_, bool fixed_step_time_) {
    fixed_step_stock = fixed_step_stock_;
    fixed_step_time = fixed_step_time_;
    if (fixed_step_stock == 1 && fixed_step_time == 1){
        fixed_dt = tn / numdiff_t;
        fixed_ds = s_max / numdiff_s;
    }

    // Code of the step size is not fixed

}


void Option::compute_solution_grid(std::string method) {
    if (type_option == "Call"  || type_option == "Put" && european == 1){
        solution_grid = solution_grid_vanilla(method);
    }


}

void Option::create_boundary_condition_time(Eigen::VectorXd& boundary_vector) {

    if (fixed_step_time == 1 && fixed_step_stock == 1){
        for (int row=0; row < boundary_vector.rows(); ++row) {
            if (type_option == "Call")
                boundary_vector(row, 0) = std::max((row + 1) * fixed_ds - strike, 0.0);
            else if (type_option == "Put")
                boundary_vector(row, 0) = std::max(strike - (row + 1) * fixed_ds, 0.0);
        }
    }

}


Eigen::MatrixXd Option::solution_grid_vanilla(std::string method) {
    Eigen::MatrixXd solution_matrix;

    if (method == "Implicit"){
        // Implementation of the solution grid computation with Dirichlet boundary conditions
        if (stock_boundary_condition == "Dirichlet"){

            Eigen::MatrixXd implicit_matrix = implicit_matrix_vanilla().inverse();
            Eigen::VectorXd current_vt;
            Eigen::VectorXd previous_vt;

            Eigen::Index implicit_matrix_rows = implicit_matrix.rows();

            current_vt.resize(implicit_matrix_rows);
            create_boundary_condition_time(current_vt);
            previous_vt.resize(implicit_matrix_rows);

            Eigen::Index solution_matrix_cols = numdiff_t + 1;

            solution_matrix.resize(implicit_matrix_rows, solution_matrix_cols);

            Eigen::VectorXd adjustment_dirichlet = Eigen::VectorXd::Zero(implicit_matrix_rows);

            double a = parameter_a(1);
            double c = parameter_c(numdiff_s - 1);

            for (int t = 0; t < solution_matrix_cols; ++t) {
                // Add adjustment for Dirichlet
                if (type_option == "Call"){
                    adjustment_dirichlet(implicit_matrix_rows - 1,0) = c * (s_max - strike * exp( - deterministic_rate * (t - 1) * fixed_dt));
                    adjustment_dirichlet(0,0) = 0;
                } else if (type_option == "Put"){
                    adjustment_dirichlet(implicit_matrix_rows - 1,0) = 0;
                    adjustment_dirichlet(0,0) = a * strike * exp( - deterministic_rate * (t - 1) * fixed_dt);
                }

                // Perform the iteration
                if (t == 0) {
                    solution_matrix.col(solution_matrix_cols - 1) = current_vt;
                } else {
                    previous_vt = implicit_matrix * (current_vt - adjustment_dirichlet);
                    solution_matrix.col(solution_matrix_cols - t - 1) = previous_vt;
                    current_vt = previous_vt;
                }
            }

            Eigen::MatrixXd solution_matrix_adjusted;
            solution_matrix_adjusted.resize(implicit_matrix_rows + 2, solution_matrix_cols);
            Eigen::MatrixXd & ref_solution_matrix = solution_matrix;
            solution_matrix_adjusted.middleRows(1,implicit_matrix_rows) = ref_solution_matrix;

            for (int col = 0; col < solution_matrix_cols; ++col) {
                if (type_option == "Call"){
                    solution_matrix_adjusted(0,solution_matrix_cols - col - 1) = 0;
                    solution_matrix_adjusted(implicit_matrix_rows + 1, solution_matrix_cols - col - 1) = s_max - strike * exp( - deterministic_rate * col * fixed_dt);
                } else if (type_option == "Put") {
                    solution_matrix_adjusted(0,solution_matrix_cols - col - 1) = strike * exp( - deterministic_rate * col * fixed_dt);
                    solution_matrix_adjusted(implicit_matrix_rows + 1, solution_matrix_cols - col - 1) = 0;
                }
            }
            return solution_matrix_adjusted;
        }

        if (stock_boundary_condition == "Neumann"){

            double a = parameter_a(1);
            double c = parameter_c(numdiff_s - 1);

            Eigen::MatrixXd implicit_matrix = implicit_matrix_vanilla();
            Eigen::Index implicit_matrix_rows = implicit_matrix.rows();
            implicit_matrix(0,0) += 2 * a;
            implicit_matrix(0,1) -= a;

            implicit_matrix(implicit_matrix_rows - 1, implicit_matrix_rows - 1) += 2 * c;
            implicit_matrix(implicit_matrix_rows - 1, implicit_matrix_rows - 2) -= c;

            implicit_matrix = implicit_matrix.inverse();

            Eigen::VectorXd current_vt;
            current_vt.resize(implicit_matrix_rows);
            create_boundary_condition_time(current_vt);

            Eigen::VectorXd previous_vt;
            previous_vt.resize(implicit_matrix_rows);

            Eigen::Index solution_matrix_cols = numdiff_t + 1;

            solution_matrix.resize(implicit_matrix_rows, solution_matrix_cols);

            for (int t = 0; t < solution_matrix_cols; ++t) {
                // Perform the iteration
                if (t == 0) {
                    solution_matrix.col(solution_matrix_cols - 1) = current_vt;
                } else {
                    previous_vt = implicit_matrix * current_vt;
                    solution_matrix.col(solution_matrix_cols - t - 1) = previous_vt;
                    current_vt = previous_vt;
                }
            }

            Eigen::MatrixXd solution_matrix_adjusted;
            solution_matrix_adjusted.resize(implicit_matrix_rows + 2, solution_matrix_cols);
            Eigen::MatrixXd & ref_solution_matrix = solution_matrix;
            solution_matrix_adjusted.middleRows(1,implicit_matrix_rows) = ref_solution_matrix;

            for (int col = 0; col < solution_matrix_cols; ++col) {
                Eigen::Index column_number = solution_matrix_cols - col - 1;
                Eigen::Index last_row = implicit_matrix_rows + 1;
                solution_matrix_adjusted(0,column_number) = 2 * solution_matrix_adjusted(1,column_number)
                            - solution_matrix_adjusted(2,column_number);
                solution_matrix_adjusted(last_row, column_number) = 2 * solution_matrix_adjusted(last_row - 1, column_number)
                            - solution_matrix_adjusted(last_row - 2, column_number);

            }
            return solution_matrix_adjusted;
        }}

    if (method == "Crank-Nicholson"){
        if (stock_boundary_condition == "Dirichlet"){

            Eigen::MatrixXd implicit_matrix = implicit_matrix_vanilla();
            Eigen::MatrixXd implicit_matrix_cn = - implicit_matrix;
            Eigen::Index implicit_matrix_rows = implicit_matrix_cn.cols();
            Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(implicit_matrix_rows, implicit_matrix_rows);
            const Eigen::MatrixXd & ref_identity = identity;
            implicit_matrix_cn += 3 * ref_identity;
            implicit_matrix += ref_identity;
            implicit_matrix = implicit_matrix.inverse();


            Eigen::VectorXd current_vt;
            Eigen::VectorXd previous_vt;

            current_vt.resize(implicit_matrix_rows);
            create_boundary_condition_time(current_vt);
            previous_vt.resize(implicit_matrix_rows);

            Eigen::Index solution_matrix_cols = numdiff_t + 1;
            solution_matrix.resize(implicit_matrix_rows, solution_matrix_cols);

            Eigen::VectorXd adjustment_dirichlet = Eigen::VectorXd::Zero(implicit_matrix_rows);

            double a = parameter_a(1);
            double c = parameter_c(numdiff_s - 1);

            for (int t = 0; t < solution_matrix_cols; ++t) {
                // Add adjustment for Dirichlet
                if (type_option == "Call"){
                    adjustment_dirichlet(adjustment_dirichlet.rows() - 1,0) = c * (2 * s_max - strike * exp( - deterministic_rate * (t - 1) * fixed_dt)
                            - strike * exp( - deterministic_rate * std::max(t - 2,0) * fixed_dt));
                    adjustment_dirichlet(0,0) = 0;
                } else if (type_option == "Put"){
                    adjustment_dirichlet(adjustment_dirichlet.rows() - 1,0) = 0;
                    adjustment_dirichlet(0,0) = a * strike * ( exp( - deterministic_rate * (t - 1) * fixed_dt)
                            + exp( - deterministic_rate * std::max(t - 2,0) * fixed_dt));
                }

                // Perform the iteration
                if (t == 0) {
                    solution_matrix.col(solution_matrix_cols - 1) = current_vt;
                } else {
                    previous_vt = implicit_matrix * (implicit_matrix_cn * current_vt - adjustment_dirichlet);
                    solution_matrix.col(solution_matrix_cols - t - 1) = previous_vt;
                    current_vt = previous_vt;
                }
            }

            Eigen::MatrixXd solution_matrix_adjusted;
            solution_matrix_adjusted.resize(implicit_matrix_rows + 2, solution_matrix_cols);
            Eigen::MatrixXd & ref_solution_matrix = solution_matrix;
            solution_matrix_adjusted.middleRows(1,implicit_matrix_rows) = ref_solution_matrix;

            for (int col = 0; col < solution_matrix_cols; ++col) {
                if (type_option == "Call"){
                    solution_matrix_adjusted(0,solution_matrix_cols - col - 1) = 0;
                    solution_matrix_adjusted(implicit_matrix_rows + 1, solution_matrix_cols - col - 1) = s_max - strike * exp( - deterministic_rate * col * fixed_dt);
                } else if (type_option == "Put") {
                    solution_matrix_adjusted(0,solution_matrix_cols - col - 1) = strike * exp( - deterministic_rate * col * fixed_dt);
                    solution_matrix_adjusted(implicit_matrix_rows + 1, solution_matrix_cols - col - 1) = 0;
                }
            }
            return solution_matrix_adjusted;
        }
        if (stock_boundary_condition == "Neumann"){
            Eigen::MatrixXd implicit_matrix = implicit_matrix_vanilla();
            Eigen::MatrixXd implicit_matrix_cn = - implicit_matrix;
            Eigen::Index implicit_matrix_rows = implicit_matrix_cn.rows();
            Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(implicit_matrix_rows, implicit_matrix_rows);
            const Eigen::MatrixXd & ref_identity = identity;
            implicit_matrix_cn += 3 * ref_identity;
            implicit_matrix += ref_identity;

            double a = parameter_a(1);
            double c = parameter_c(numdiff_s - 1);

            implicit_matrix(0,0) += 2 * a;
            implicit_matrix(0,1) -= a;

            implicit_matrix(implicit_matrix.rows() - 1, implicit_matrix.cols() - 1) += 2 * c;
            implicit_matrix(implicit_matrix.rows() - 1, implicit_matrix.cols() - 2) -= c;

            implicit_matrix_cn(0,0) -= 2 * a;
            implicit_matrix_cn(0,1) += a;

            implicit_matrix_cn(implicit_matrix.rows() - 1, implicit_matrix.cols() - 1) -= 2 * c;
            implicit_matrix_cn(implicit_matrix.rows() - 1, implicit_matrix.cols() - 2) += c;

            implicit_matrix = implicit_matrix.inverse();

            Eigen::VectorXd current_vt;
            current_vt.resize(implicit_matrix_rows);
            create_boundary_condition_time(current_vt);

            Eigen::VectorXd previous_vt;
            previous_vt.resize(implicit_matrix_rows);

            Eigen::Index solution_matrix_cols = numdiff_t + 1;
            solution_matrix.resize(implicit_matrix_rows, solution_matrix_cols);

            for (int t = 0; t < solution_matrix_cols; ++t) {
                // Perform the iteration
                if (t == 0) {
                    solution_matrix.col(solution_matrix_cols - 1) = current_vt;
                } else {
                    previous_vt = implicit_matrix * (implicit_matrix_cn * current_vt);
                    solution_matrix.col(solution_matrix_cols - t - 1) = previous_vt;
                    current_vt = previous_vt;
                }
            }

            Eigen::MatrixXd solution_matrix_adjusted;

            solution_matrix_adjusted.resize(implicit_matrix_rows + 2, solution_matrix_cols);
            Eigen::MatrixXd & ref_solution_matrix = solution_matrix;
            solution_matrix_adjusted.middleRows(1,implicit_matrix_rows) = ref_solution_matrix;

            for (int col = 0; col < solution_matrix_cols; ++col) {
                Eigen::Index column_number = solution_matrix_cols - col - 1;
                Eigen::Index last_row = implicit_matrix_rows + 1;
                solution_matrix_adjusted(0,column_number) = 2 * solution_matrix_adjusted(1,column_number)
                                                            - solution_matrix_adjusted(2,column_number);
                solution_matrix_adjusted(last_row, column_number) = 2 * solution_matrix_adjusted(last_row - 1, column_number)
                                                                    - solution_matrix_adjusted(last_row - 2, column_number);
            }
            return solution_matrix_adjusted;
        }

    }

}


void Option::print_solution_grid(int precision) {
    std::cout << "Type of Option: " + type_option + "\n";
    std::cout << "Stock boundary condition: " + stock_boundary_condition + "\n\n";
    std::cout << "Legend of the matrix: s/t 1 2 . . T\n";
    std::cout << "                      0   x x x x x\n";
    std::cout << "                      1   x x x x x\n";
    std::cout << "                      .   x x x x x\n";
    std::cout << "                      .   x x x x x\n";
    std::cout << "                  S_max   x x x x x\n\n";
    std::cout << "Printing the solution grid:\n\n";

    Eigen::MatrixXd printing_matrix;
    printing_matrix.resize(solution_grid.rows() + 1, solution_grid.cols() + 1);
    printing_matrix.bottomRightCorner(solution_grid.rows(),solution_grid.cols()) = solution_grid;

    if (fixed_step_time == 1){
        for (int col = 1; col < printing_matrix.cols(); ++col) {
            printing_matrix(0, col) = fixed_dt * (col - 1);
        }

    }

    if (fixed_step_stock == 1){
        for (int row = 1; row < printing_matrix.rows(); ++row) {
            printing_matrix(row, 0) = fixed_ds * (row - 1);
        }
    }

    std::streamsize prec = std::cout.precision();
    std::cout << std::setprecision(precision) << printing_matrix << std::setprecision(prec) << "\n";
}
