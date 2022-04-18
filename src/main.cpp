#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include "Option.h"
#include <cmath>
#include "Regression.h"
#include "FileManagement.h"
#include <chrono>
#include "nonuniform_grid.h"
#include <random>
#include "Simulation.h"
#include <tuple>
#include "HestonOption.h"
#include "OptimalExecution.h"




int main() {

  // Generation of option prices
  // auto o = Option("Call");
  // o.set_tn(1);
  // o.set_s_max(100);
  // o.set_numdiff_t(10);
  // o.set_numdiff_s(10);
  // o.set_volatility(.1);
  // o.set_interest_rate(.15);
  // o.set_strike(50);
  // o.set_stock_boundary_condition("Dirichlet");
  // o.fixed_difference_step();
  //
  // o.compute_solution_grid("Crank-Nicholson");
  // o.print_solution_grid();


  // Heston option Monte-Carlo Simulation
  // auto h = HestonOption("Call", 1, 10, 100, 50);
  // h.set_params_stock_process(45, .1);
  // h.set_prams_variance_process(.2, 1, 1, .25, .5);
  // auto pr = h.compute_price();
  //
  // std::cout << pr << "\n";


  // Generation of optimal execution paths
   // double a_coeff_, double b_coeff_, double sigma_coeff_, double k_coeff_, double phi_coeff_
   // auto o = OptimalExecution(.6,.1, 0.1,.1,.6);
   // // num_diff_t_, int num_diff_q_
   // o.set_num_steps(100,100);
   // // true_impact, trading_speed, impact_nonlinearity
   // o.compute_liquidation("Linear", 1);
   // Eigen::VectorXd h = o.get_value_process();
   // std::cout << "Value Process\n" << h << "\n";
   // Eigen::VectorXd j = o.get_optimal_speed_process();
   // std::cout << "Speed Process\n" << j << "\n";
   // Eigen::VectorXd s = o.get_stock_process();
   // std::cout << "Stock Process\n" << s << "\n";
   // Eigen::VectorXd c = o.get_cash_process();
   // std::cout << "Cash Process\n" << c << "\n";
   // Eigen::VectorXd i = o.get_inventory_process();
   // std::cout << "Inventory Process\n" << i << "\n";
   //
   // std::cout << o.get_value_matrix();

    // Generate 500 optimal execution scenarios and save it in the csv files
    // Eigen::MatrixXd nonlinear_speed_csv_value,nonlinear_speed_csv_optimal_speed, nonlinear_speed_csv_stock, nonlinear_speed_csv_cash, nonlinear_speed_csv_inventory;
    // nonlinear_speed_csv_value.resize(101,100); nonlinear_speed_csv_optimal_speed.resize(101,100); nonlinear_speed_csv_stock.resize(101,100); nonlinear_speed_csv_cash.resize(101,100); nonlinear_speed_csv_inventory.resize(101,100);
    // Eigen::MatrixXd linear_speed_csv_value,linear_speed_csv_optimal_speed, linear_speed_csv_stock, linear_speed_csv_cash, linear_speed_csv_inventory;
    // linear_speed_csv_value.resize(101,100); linear_speed_csv_optimal_speed.resize(101,100); linear_speed_csv_stock.resize(101,100); linear_speed_csv_cash.resize(101,100); linear_speed_csv_inventory.resize(101,100);
    // for (int sim = 0; sim < 100; ++sim) {
        // Trading speed: Nonlinear, true impact: Nonlinear
        // auto o = OptimalExecution(.08,.06, .1,.08,.06);
        // o.set_num_steps(100,100);
        // o.compute_liquidation("Nonlinear", "Nonlinear", .6);
        // Eigen::VectorXd value_process = o.get_value_process();
        // nonlinear_speed_csv_value.col(sim) = value_process;
        // Eigen::VectorXd optimal_speed_process = o.get_optimal_speed_process();
        // nonlinear_speed_csv_optimal_speed.col(sim) = optimal_speed_process;
        // Eigen::VectorXd stock_process = o.get_stock_process();
        // nonlinear_speed_csv_stock.col(sim) = stock_process;
        // Eigen::VectorXd cash_process = o.get_cash_process();
        // nonlinear_speed_csv_cash.col(sim) = cash_process;
        // Eigen::VectorXd inventory_process = o.get_inventory_process();
        // nonlinear_speed_csv_inventory.col(sim) = inventory_process;

        // Trading speed: Linear, true impact: Nonlinear
        // auto ol = OptimalExecution(.6,.1, 0.1,.1,.6);
        // ol.set_num_steps(100,100);
        // ol.compute_liquidation("Linear", 1);
        // Eigen::VectorXd value_process_lin = ol.get_value_process();
        // linear_speed_csv_value.col(sim) = value_process_lin;
        // Eigen::VectorXd optimal_speed_process_lin = ol.get_optimal_speed_process();
        // linear_speed_csv_optimal_speed.col(sim) = optimal_speed_process_lin;
        // Eigen::VectorXd stock_process_lin = ol.get_stock_process();
        // linear_speed_csv_stock.col(sim) = stock_process_lin;
        // Eigen::VectorXd cash_process_lin = ol.get_cash_process();
        // linear_speed_csv_cash.col(sim) = cash_process_lin;
        // Eigen::VectorXd inventory_process_lin = ol.get_inventory_process();
        // linear_speed_csv_inventory.col(sim) = inventory_process_lin;
    // }
    // saveData("data/nonlinear_speed_csv_value.csv", nonlinear_speed_csv_value);
    // saveData("data/nonlinear_speed_csv_optimal_speed.csv", nonlinear_speed_csv_optimal_speed);
    // saveData("data/nonlinear_speed_csv_stock.csv", nonlinear_speed_csv_stock);
    // saveData("data/nonlinear_speed_csv_cash.csv", nonlinear_speed_csv_cash);
    // saveData("data/nonlinear_speed_csv_inventory.csv", nonlinear_speed_csv_inventory);

    // saveData("data/linear_speed_csv_value.csv", linear_speed_csv_value);
    // saveData("data/linear_speed_csv_optimal_speed.csv", linear_speed_csv_optimal_speed);
    // saveData("data/linear_speed_csv_stock.csv", linear_speed_csv_stock);
    // saveData("data/linear_speed_csv_cash.csv", linear_speed_csv_cash);
    // saveData("data/linear_speed_csv_inventory.csv", linear_speed_csv_inventory);





    return 0;
}
