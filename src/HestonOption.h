//
// Created by Vlad Argunov on 10/01/2022.
//

#ifndef QUANTKIT_HESTONOPTION_H
#define QUANTKIT_HESTONOPTION_H

#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "Simulation.h"


class HestonOption {
private:
    Simulation sim;

    std::string type_option; // "Call", "Put"
    double tn;
    int n_processes;
    int n_steps;

    double stock_price_start;
    double strike;

    double mu;
    double variance_start;
    double kappa;
    double theta;
    double sigma;
    double rho;

public:
    explicit HestonOption(std::string type_option_, double tn_, int n_processes_, int n_steps_, double strike_);
    void set_params_stock_process(double stock_price_start_, double mu_);
    void set_prams_variance_process(double variance_start_, double kappa_, double theta_, double sigma_, double rho_);
    double compute_price();
};


#endif //QUANTKIT_HESTONOPTION_H
