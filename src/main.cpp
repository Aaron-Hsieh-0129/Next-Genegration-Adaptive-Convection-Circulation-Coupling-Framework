#include "../2DVVM/src/Declare.hpp"
#include "../CSSWM/src/construction.hpp"
#include "../CSSWM/src/define.hpp"
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <netcdf>

using namespace netCDF;

// #define AB2_COUPLE


// CASE0: Nothing, CASE1:Bubble
Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double Kx, double Kz) {
    return Config_VVM(3.0, 200.0, 200.0, 100000, 20000, 1500000.0, 
                      10000, path, 10, 
                      Kx, Kz, 0.01, 1E-22, 9.80665, 1003.5, 716.5, 287.0, 
                      2.5E6, 1E5, 96500.0, addforcingtime, CASE);
}



// void output_qall(std::string dir,int n, double q[6][model_csswm.nx][model_csswm.ny]);
// void output_Qall(std::string dir,int n, double Q[6][model_csswm.nx][model_csswm.ny]);

vvm**** allocate_and_initialize(int dim1, int dim2, int dim3) {
    // Allocate memory for 3D array (layers x model_csswm.nx x model_csswm.ny)
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
    // Allocate memory for 3D array (layers x model_csswm.nx x model_csswm.ny)
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
}

void deallocate(vvm**** array, int dim1, int dim2, int dim3) {
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
}

struct vvm_index {
    int p, i, j;
};


int main(int argc, char **argv) {
    omp_set_num_threads(128);
    Eigen::setNbThreads(1);

    std::string path = "/data/Aaron/TMIF_AB2/newTur_dt300_1_cloud__csswm_2E5diff_p1/";

    Config_CSSWM config_csswm(300., 1., 1., 0.1, 86400 * 3 * 24., path + "csswm/", 
                        1, 2E5, 2E5, 0.06, 1200. * 60.);
    CSSWM model_csswm(config_csswm);

    CSSWM::Init::Init2d(model_csswm);
    

    Config_VVM**** config_vvms = allocate_and_initialize_config(6, model_csswm.nx, model_csswm.ny);
    std::string path_vvm;

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < model_csswm.nx; i++) {
            for (int j = 0; j < model_csswm.ny; j++) {
                path_vvm = path + "vvm/" + std::to_string(p) + "_" + std::to_string(i) + "_" + std::to_string(j) + "/";
                // if (p == 1 && (i >= model_csswm.nx/2-5 && i <= model_csswm.nx/2+5) && j == model_csswm.ny/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (i >= model_csswm.nx/2-2 && i <= model_csswm.nx/2+2) && j == model_csswm.ny/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && i == model_csswm.nx/2 && j == model_csswm.ny/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (i >= model_csswm.nx/2-4 && i <= model_csswm.nx/2+4) && (j >= model_csswm.nx/2-4 && j <= model_csswm.nx/2+4)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (model_csswm.nx/2-3 <= i && i <= model_csswm.nx/2+3) && (j == model_csswm.ny/2)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                if (p == 1 && (model_csswm.ny/2-10 <= j && j <= model_csswm.ny/2+10) && (i == model_csswm.nx/2)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                else config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 0, 70, 70));
            }
        }
    }
    printf("Configurations are set.\n");

    int total_size = 81;
    vvm_index vvms_index[total_size];
    int count = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i <= model_csswm.nx-2; i++) {
            for (int j = 2; j <= model_csswm.ny-2; j++) {
                if (p == 1 && (model_csswm.nx/2-30 <= i && i <= model_csswm.nx/2+30) && j == model_csswm.ny/2) {
                    vvms_index[count] = {p, i, j};
                    count++;
                }
                if (p == 1 && i == model_csswm.nx/2 && (model_csswm.ny/2-10 <= j && j <= model_csswm.ny/2+10 && j != model_csswm.ny/2)) {
                    vvms_index[count] = {p, i, j};
                    count++;
                }

                // if (p == 1 && (i >= model_csswm.nx/2-10 && i <= model_csswm.nx/2+10) && (j >= model_csswm.nx/2-10 && j <= model_csswm.nx/2+10)) {
                //     vvms_index[count] = {p, i, j};
                //     count++;
                // }
            }
        }
    }
    printf("count: %d\n", count);
    if (count != total_size) {
        printf("Error: count != total_size\n");
        return 1;
    }
    
    vvm**** vvms = allocate_and_initialize(6, model_csswm.nx, model_csswm.ny);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;
        vvms[p][i][j] = new vvm(*config_vvms[p][i][j]);
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif
    printf("VVMs are declared.\n");

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;
        #if defined(LOADFROMPREVIOUSFILE)
            vvm::Init::LoadFromPreviousFile(*vvms[p][i][j]);
        #elif defined(LOAD2DINIT)
            vvm::Init::Load2DInit(*vvms[p][i][j]);
        #else
            vvm::Init::Init1d(*vvms[p][i][j]);
            vvm::Init::Init2d(*vvms[p][i][j]);
            #ifndef PETSC
                vvm::PoissonSolver::InitPoissonMatrix(*vvms[p][i][j]);
            #endif
        #endif
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif
    printf("VVMs are initialized.\n");

    int n_csswm = 0;
    double temp_csswm = model_csswm.timeend / model_csswm.dt, temp_vvm = model_csswm.timeend / config_vvms[1][model_csswm.nx/2][model_csswm.ny/2]->dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    CSSWM::Outputs::create_all_directory(model_csswm.outputpath);
    // create Q_all directory
    CSSWM::Outputs::create_directory(model_csswm.outputpath + "Q_all/");

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;
        vvm::Output::create_all_directory(*vvms[p][i][j]);
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif
    vvm::Output::create_directory(path + "vvm/q_all/");

    #ifdef NCOUTPUT
        CSSWM::Outputs::grid_nc(model_csswm);
    #endif
    #ifdef TXTOUTPUT
        CSSWM::Outputs::grid(model_csswm);
    #endif

    int vvm_nx = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->nx;
    int vvm_nz = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->nz;

    double exchange_coeff = 287. / 9.80665;
    double Q = 0.;
    
    double coupling_csswm_param = 1.;
    double coupling_vvm_param = 1.;

    double thm_mean = 0.;
    double th_mean = 0.;
    double th_mean_all[6][model_csswm.nx][model_csswm.ny];
    double Q_all[6][model_csswm.nx][model_csswm.ny];
    #if defined(AB2_COUPLE)
        double q_all[2][6][model_csswm.nx][model_csswm.ny];
    #else
        double q_all[6][model_csswm.nx][model_csswm.ny];
    #endif
    // initialize Q_all, q_all
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < model_csswm.nx; i++) {
            for (int j = 0; j < model_csswm.ny; j++) {
                Q_all[p][i][j] = 0.;
                th_mean_all[p][i][j] = 0.;
                #if defined(AB2_COUPLE)
                    q_all[0][p][i][j] = 0.;
                    q_all[1][p][i][j] = 0.;
                #else
                    q_all[p][i][j] = 0.;
                #endif
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif

    while (vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step < nmax_vvm || n_csswm < nmax_csswm) {

        double time_vvm = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step * vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->dt;
        double time_csswm = n_csswm * model_csswm.dt;
        printf("n_vvm: %d, time_vvm: %f, n_csswm: %d,  time_csswm: %f\n", vvms[1][47][47]->step, time_vvm, n_csswm, time_csswm);

        if (time_vvm == time_csswm) {
            if (n_csswm % model_csswm.outputstep == 0 || n_csswm == model_csswm.timeend-1 || n_csswm == model_csswm.timeend-2) {
                #ifdef NCOUTPUT
                    CSSWM::Outputs::huv_nc(n_csswm, model_csswm);
                #endif
            }


            CSSWM::Iteration::ph_pt_4(model_csswm);
            #ifndef Advection
                CSSWM::Iteration::pu_pt_4(model_csswm);
                CSSWM::Iteration::pv_pt_4(model_csswm);
            #endif

            // Boundary exchange and interpolation
            model_csswm.BP_h(model_csswm);
            #ifndef Advection
                model_csswm.BP_wind_interpolation2(model_csswm);
            #endif

            #if defined(DIFFUSION)
                CSSWM::NumericalProcess::DiffusionAll(model_csswm);
            #endif

            model_csswm.BP_h(model_csswm);
            model_csswm.BP_wind_interpolation2(model_csswm);

            
            #if defined(TIMEFILTER) && !defined(AB2)
                CSSWM::NumericalProcess::timeFilterAll(model_csswm);
            #endif

            // next step
            // CSSWM::Iteration::nextTimeStep(model_csswm);
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
     
        if (vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step % vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->OUTPUTSTEP == 0) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int size = 0; size < total_size; size++) {
                int p = vvms_index[size].p;
                int i = vvms_index[size].i;
                int j = vvms_index[size].j;
                #if defined(OUTPUTTXT)
                    vvm::Output::outputalltxt(vvms[p][i][j]->step, *vvms[p][i][j]);
                #endif

                #if defined(OUTPUTNC)
                    vvm::Output::output_nc(vvms[p][i][j]->step, *vvms[p][i][j]);
                #endif
            }
            #ifdef _OPENMP
            #pragma omp barrier
            #endif

            // #if defined(AB2_COUPLE)
            //     output_qall(path + "vvm/q_all/", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, q_all[1]);
            // #else
            //     output_qall(path + "vvm/q_all/", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, q_all);
            // #endif
        }

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;

            vvm::Iteration::pzeta_pt(*vvms[p][i][j]);
            vvm::Iteration::pth_pt(*vvms[p][i][j]);
            #if defined(WATER)
                vvm::Iteration::pqv_pt(*vvms[p][i][j]);
                vvm::Iteration::pqc_pt(*vvms[p][i][j]);
                vvm::Iteration::pqr_pt(*vvms[p][i][j]);

                if (vvms[p][i][j]->step * vvms[p][i][j]->dt <= vvms[p][i][j]->addforcingtime) vvms[p][i][j]->status_for_adding_forcing = true;
                else vvms[p][i][j]->status_for_adding_forcing = false;

                // Generate new random th perturbation for tropical forcing case
                if (vvms[p][i][j]->status_for_adding_forcing == true) {
                    vvm::Init::RandomPerturbation(*vvms[p][i][j], vvms[p][i][j]->step);
                }
                vvm::AddForcing(*vvms[p][i][j]);
            #endif
            vvm::BoundaryProcess2D_all(*vvms[p][i][j]);

            vvm::PoissonSolver::pubarTop_pt(*vvms[p][i][j]);
            vvm::PoissonSolver::cal_w(*vvms[p][i][j], p, i, j);
            vvm::PoissonSolver::cal_u(*vvms[p][i][j]);
            
            vvm::Iteration::updateMean(*vvms[p][i][j]);
            vvm::Turbulence::RKM_RKH(*vvms[p][i][j]);

            #if defined(WATER)
                vvm::MicroPhysics::autoconversion(*vvms[p][i][j]);
                vvm::MicroPhysics::accretion(*vvms[p][i][j]);
                vvm::MicroPhysics::evaporation(*vvms[p][i][j]);
                vvm::MicroPhysics::condensation(*vvms[p][i][j]); // saturation adjustment

                // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
                vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qvp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
                vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qcp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
                vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qrp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
            #endif

            vvm::BoundaryProcess2D_all(*vvms[p][i][j]);

            #if defined(TIMEFILTER) && !defined(AB2)
                vvm::NumericalProcess::timeFilterAll(*vvms[p][i][j]);
            #endif
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        // Exchange information here
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:th_mean)
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;
            if (time_csswm == time_vvm) {
                th_mean = 0.;
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        #if defined(AB2)
                            th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                        #else
                            th_mean += vvms[p][i][j]->thm[i_vvm][k_vvm];
                        #endif
                    }
                }
                th_mean /= ((vvm_nx-2) * (vvm_nz-2));
                th_mean_all[p][i][j] = th_mean;
            }
            #if defined(AB2_COUPLE)
                q_all[(vvms[p][i][j]->step+1)%2][p][i][j] = coupling_vvm_param * (model_csswm.hp[p][i][j] / exchange_coeff - th_mean) / model_csswm.dt;
                if (vvms[p][i][j]->step == 1) q_all[1][p][i][j] = q_all[0][p][i][j];
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        vvms[p][i][j]->thp[i_vvm][k_vvm] += 1.5*vvms[p][i][j]->dt*q_all[(vvms[p][i][j]->step+1)%2][p][i][j] - 0.5*vvms[p][i][j]->dt*q_all[vvms[p][i][j]->step%2][p][i][j];
                    }
                }
            #else
                q_all[p][i][j] = coupling_vvm_param * (model_csswm.hp[p][i][j] / exchange_coeff - th_mean_all[p][i][j]) / model_csswm.dt;
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        vvms[p][i][j]->thp[i_vvm][k_vvm] += vvms[p][i][j]->dt * q_all[p][i][j];
                    }
                }
            #endif
            
            
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif


        if (time_csswm == time_vvm) {
            // Exchange information here
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int size = 0; size < total_size; size++) {
                int p = vvms_index[size].p;
                int i = vvms_index[size].i;
                int j = vvms_index[size].j;
                
                Q_all[p][i][j] = (exchange_coeff * th_mean_all[p][i][j] - model_csswm.h[p][i][j]) / model_csswm.dt;

                model_csswm.hp[p][i][j] += coupling_csswm_param * Q_all[p][i][j] * model_csswm.dt;
            }
            #ifdef _OPENMP
            #pragma omp barrier
            #endif

            // output_Qall(model_csswm.outputpath + (std::string) "Q_all/", n_csswm, Q_all);

            CSSWM::Iteration::nextTimeStep(model_csswm);

            n_csswm++;
        }

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;
            vvm::Iteration::nextTimeStep(*vvms[p][i][j]);
            vvms[p][i][j]->step++;
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

    }

    deallocate_config(config_vvms, 6, model_csswm.nx, model_csswm.ny);
    deallocate(vvms, 6, model_csswm.nx, model_csswm.ny);
    return 0;
}

// void output_qall(std::string dir, int n, double q[6][model_csswm.nx][model_csswm.ny]) {
//     NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", model_csswm.nx);
//     NcDim yDim = dataFile.addDim("y", model_csswm.ny);
//     NcDim lonDim = dataFile.addDim("lon", model_csswm.nx);
//     NcDim latDim = dataFile.addDim("lat", model_csswm.ny);

//     std::vector<NcDim> xyDim, lonlatDim;
//     xyDim.push_back(p);
//     xyDim.push_back(xDim);
//     xyDim.push_back(yDim);

//     NcVar q_all = dataFile.addVar("q", ncDouble, xyDim);

//     std::vector<size_t> startp, countp;
//     startp.push_back(0);
//     startp.push_back(0);
//     startp.push_back(0);
//     countp.push_back(1);
//     countp.push_back(model_csswm.nx);
//     countp.push_back(model_csswm.ny);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         q_all.putVar(startp, countp, q[p]);
//     }
//     return;
// }

// void output_Qall(std::string dir,int n, double Q[6][model_csswm.nx][model_csswm.ny]) {
//     NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);       
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", model_csswm.nx);
//     NcDim yDim = dataFile.addDim("y", model_csswm.ny);
//     NcDim lonDim = dataFile.addDim("lon", model_csswm.nx);
//     NcDim latDim = dataFile.addDim("lat", model_csswm.ny);

//     std::vector<NcDim> xyDim, lonlatDim;
//     xyDim.push_back(p);
//     xyDim.push_back(xDim);
//     xyDim.push_back(yDim);

//     NcVar q_all = dataFile.addVar("q", ncDouble, xyDim);

//     std::vector<size_t> startp, countp;
//     startp.push_back(0);
//     startp.push_back(0);
//     startp.push_back(0);
//     countp.push_back(1);
//     countp.push_back(model_csswm.nx);
//     countp.push_back(model_csswm.ny);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         q_all.putVar(startp, countp, Q[p]);
//     }
//     return;
// }
