//
// Created by Vlad Argunov on 19/10/2021.
//

#ifndef QUANTKIT_FILEMANAGEMENT_H
#define QUANTKIT_FILEMANAGEMENT_H

#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <stdlib.h>
#include <vector>


// Works only with numerical data, no headers
// Ensure that the number of columns (commas) is equal in each row
// Should work fine with Excel files converted into CSV
void read_csv_into_matrix(const std::string& path, Eigen::MatrixXd & result );
std::vector<int> size_of_csv(const std::string& path);
void saveData(std::string fileName, Eigen::MatrixXd  matrix);


#endif //QUANTKIT_FILEMANAGEMENT_H
