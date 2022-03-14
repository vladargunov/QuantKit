//
// Created by Vlad Argunov on 10/01/2022.
//

#include "HestonOption.h"

#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "Simulation.h"

HestonOption::HestonOption(std::string type_option_, double tn_, int n_processes_, int n_steps_, double strike_) {
    sim = Simulation();
    type_option = type_option_;
    tn = tn_;
    n_processes = n_processes_;
    n_steps = n_steps_;
    strike = strike_;
}

void HestonOption::set_params_stock_process(double stock_price_start_, double mu_) {
    stock_price_start = stock_price_start_;
    mu = mu_;
}

void HestonOption::set_prams_variance_process(double variance_start_, double kappa_, double theta_, double sigma_,
                                              double rho_) {
    variance_start = variance_start_;
    kappa = kappa_;
    theta = theta_;
    sigma = sigma_;
    rho = rho_;
}

double HestonOption::compute_price() {
    auto [stocks, _] = sim.generate_heston_process(n_processes, tn, stock_price_start, mu, variance_start, kappa, theta, sigma, rho, n_steps);
    Eigen::Index stocks_cols = stocks.cols();
    double total_payoffs = 0;
    for (int row = 0; row < stocks.rows(); ++row) {
        if (type_option == "Call"){
            total_payoffs += std::max(stocks(row,stocks_cols - 1) - strike, 0.0);
        }
        if (type_option == "Put"){
            total_payoffs += std::max(strike - stocks(row,stocks_cols), 0.0);
        }
    }
    return exp(- mu * tn) * total_payoffs / n_processes;
}
