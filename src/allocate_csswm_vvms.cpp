#include "allocate_csswm_vvms.hpp"

// CASE: 0, Nothing; CASE: 1, Bubble; CASE: 2, Bubble with shear
Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double xrange, double zrange, double dx, double dz, double dt, double timeend, double outputstep, double vvm_moisture_nudge_time) {
    return Config_VVM(dt, dx, dz, xrange, zrange, timeend, 
                      10000, path, outputstep, 
                      100., 100., 0.01, 1E-22, 9.80665, 1003.5, 716.5, 287.0, 
                      2.5E6, 1E5, 96500.0, addforcingtime, CASE, vvm_moisture_nudge_time);
}

vvm**** allocate_and_initialize(int dim1, int dim2, int dim3) {
    // Allocate memory for 3D array (layers x NX x NY)
    vvm**** array = new vvm***[dim1];
    for (int p = 0; p < dim1; ++p) {
        array[p] = new vvm**[dim2];
        for (int i = 0; i < dim2; ++i) {
            array[p][i] = new vvm*[dim3];
            for (int j = 0; j < dim3; ++j) {
                array[p][i][j] = nullptr; // Initialize to nullptr
            }
        }
    }
    return array;
}

Config_VVM**** allocate_and_initialize_config(int dim1, int dim2, int dim3) {
    // Allocate memory for 3D array (layers x NX x NY)
    Config_VVM**** array = new Config_VVM***[dim1];
    for (int p = 0; p < 6; ++p) {
        array[p] = new Config_VVM**[dim2];
        for (int i = 0; i < dim2; ++i) {
            array[p][i] = new Config_VVM*[dim3];
            for (int j = 0; j < dim3; ++j) {
                array[p][i][j] = nullptr; // Initialize to nullptr
            }
        }
    }
    return array;
}

void deallocate_config(Config_VVM**** array, int dim1, int dim2, int dim3) {
    if (array == nullptr) return;
    
    for (int p = 0; p < dim1; ++p) {
        if (array[p] != nullptr) {
            for (int i = 0; i < dim2; ++i) {
                if (array[p][i] != nullptr) {
                    for (int j = 0; j < dim3; ++j) {
                        if (array[p][i][j] != nullptr) delete array[p][i][j];
                    }
                    delete[] array[p][i];
                }
            }
            delete[] array[p];
        }
    }
    delete[] array;
    return;
}

void deallocate(vvm**** array, int dim1, int dim2, int dim3) {
    if (array == nullptr) return;

    for (int p = 0; p < dim1; ++p) {
        for (int i = 0; i < dim2; ++i) {
            for (int j = 0; j < dim3; ++j) {
                delete array[p][i][j];
            }
            delete[] array[p][i];
        }
        delete[] array[p];
    }
    delete[] array;
    return;
}
