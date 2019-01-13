#include <iostream>
#include "utility.h"
#include <vector>
#include "ArgumentHandler.h"
#include "data.h"
#include "Bandwidth.h"
#include "kde.h"
#include "cdc_statstics_method.h"

int main(int argc, char **argv) {
    try {
        // Handle command line arguments
        ArgumentHandler arg_handler(argc, argv);
        if (arg_handler.processArguments() != 0) {
            return 0;
        }
        arg_handler.checkArguments();

        Data data = Data(arg_handler.filename, arg_handler.x_index, arg_handler.y_index, arg_handler.z_index,
                         arg_handler.distance);
        Bandwidth bandwidth = Bandwidth(data.getConditional_var_num());
        std::vector<std::vector<double>> data_z = data.getData_z();
        std::vector<std::vector<double>> bandwidth_matrix = bandwidth.bandwidth_selection(data_z);

        KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(data_z, bandwidth_matrix, 1);
        kernelDensityEstimation.compute_kernel_density_estimate();
        std::vector<std::vector<double> > kernel = kernelDensityEstimation.get_kernel_density_estimate();

        std::vector<std::vector<double>> distance_x = data.getDistance_x();
        std::vector<std::vector<double>> distance_y = data.getDistance_y();
        CDCStatsticsMethod cDCStatsticsMethod = CDCStatsticsMethod();
        cDCStatsticsMethod.conduct_cdc_test(distance_x, distance_y, kernel, arg_handler.num_permutation, 1, 1);

        std::cout << "CDC statistics: " << cDCStatsticsMethod.getCdc_statistic()[0] << std::endl;
        std::cout << "p-value: " << cDCStatsticsMethod.getP_value() << std::endl;

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << " CDC will EXIT now." << std::endl;
        return -1;
    }

    return 0;
}