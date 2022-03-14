//
// Created by Vlad Argunov on 19/10/2021.
//

#include "FileManagement.h"
#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <stdlib.h>


std::vector<int> size_of_csv(const std::string& path){
// The solution is from https://stackoverflow.com/questions/415515/how-can-i-read-and-manipulate-csv-file-data-in-c?noredirect=1&lq=1
// It is enhanced to compute the rows and columns of csv
    int rows = 0;
    int cols = 0;

    std::ifstream  data(path);
    std::string line;
    bool compute_cols = 0;
    while(std::getline(data,line))
    {
        if (compute_cols == 0){
            cols = count(line.begin(), line.end(), ',');
            compute_cols = 1;
        }

        rows++;
    }
    std::vector<int> v {rows, cols + 1};
    return v;
}


void read_csv_into_matrix(const std::string& path, Eigen::MatrixXd & result ){
// The solution is from https://stackoverflow.com/questions/415515/how-can-i-read-and-manipulate-csv-file-data-in-c?noredirect=1&lq=1
// It is enhanced to handle Eigen matrices

    std::ifstream data(path);
    std::string line;

    std::vector<int> matrix_size = size_of_csv(path);
    result.resize(matrix_size[0], matrix_size[1]);

    int row = 0;
    while(std::getline(data,line))
    {

        std::stringstream  lineStream(line);
        std::string        cell;

        int col = 0;
        while(std::getline(lineStream,cell,','))
        {

            if (row == 0 && col == 0){
                cell = cell.substr(3);

            }

            double numeric_cell = atof(cell.c_str());
            result(row, col) = numeric_cell;
            col++;
        }
        row++;
    }

}

void saveData(std::string fileName, Eigen::MatrixXd  matrix)
{
    // Code from https://aleksandarhaber.com/eigen-matrix-library-c-tutorial-saving-and-loading-data-in-from-a-csv-file/
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
