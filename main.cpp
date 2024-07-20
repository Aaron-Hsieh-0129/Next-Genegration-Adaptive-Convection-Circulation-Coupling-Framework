#include "2DVVM/src/Declare.hpp"
#include "CSSWM/src/construction.hpp"
#include "CSSWM/src/define.hpp"
#include "Eigen/src/Core/products/Parallelizer.h"
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <netcdf>

using namespace netCDF;

#define AB2_COUPLE

CSSWM model_csswm;

// CASE0: Nothing, CASE1:Bubble
Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double Kx, double Kz) {
    return Config_VVM(3.0, 200.0, 200.0, 100000, 20000, 1500000.0, 
                      10000, path, 10, 
                      Kx, Kz, 0.01, 0.0, 0.0, 1E-22, 9.80665, 1003.5, 716.5, 287.0, 
                      2.5E6, 1E5, 96500.0, addforcingtime, CASE);
}

void output_qall(std::string dir,int n, double q[6][NX][NY]);
void output_Qall(std::string dir,int n, double Q[6][NX][NY]);

vvm**** allocate_and_initialize() {
    // Allocate memory for 3D array (layers x NX x NY)
    vvm**** array = new vvm***[6];
    for (int p = 0; p < 6; ++p) {
        array[p] = new vvm**[NX];
        for (int i = 0; i < NX; ++i) {
            array[p][i] = new vvm*[NY];
            for (int j = 0; j < NY; ++j) {
                array[p][i][j] = nullptr; // Initialize to nullptr
            }
        }
    }
    return array;
}

Config_VVM**** allocate_and_initialize_config() {
    // Allocate memory for 3D array (layers x NX x NY)
    Config_VVM**** array = new Config_VVM***[6];
    for (int p = 0; p < 6; ++p) {
        array[p] = new Config_VVM**[NX];
        for (int i = 0; i < NX; ++i) {
            array[p][i] = new Config_VVM*[NY];
            for (int j = 0; j < NY; ++j) {
                array[p][i][j] = nullptr; // Initialize to nullptr
            }
        }
    }
    return array;
}

void deallocate_config(Config_VVM**** array) {
    for (int p = 0; p < 6; ++p) {
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                delete array[p][i][j];
            }
            delete[] array[p][i];
        }
        delete[] array[p];
    }
    delete[] array;
}

void deallocate(vvm**** array) {
    for (int p = 0; p < 6; ++p) {
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
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

    CSSWM::Init::Init2d(model_csswm);
    
    std::string path = "/data/Aaron/TMIF/newTur_dt_4.55_cloud_2_csswm_1E6diff_60p1/";

    Config_VVM**** config_vvms = allocate_and_initialize_config();
    std::string path_vvm;

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                path_vvm = path + "vvm/" + std::to_string(p) + "_" + std::to_string(i) + "_" + std::to_string(j) + "/";
                // if (p == 1 && (i >= NX/2-5 && i <= NX/2+5) && j == NY/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (i >= NX/2-2 && i <= NX/2+2) && j == NY/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && i == NX/2 && j == NY/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (i >= NX/2-4 && i <= NX/2+4) && (j >= NX/2-4 && j <= NX/2+4)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                // if (p == 1 && (NX/2-3 <= i && i <= NX/2+3) && (j == NY/2)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                if (p == 1 && (NY/2-10 <= j && j <= NY/2+10) && (i == NX/2)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 1, 200, 200));
                else config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 0, 70, 70));
            }
        }
    }
    printf("Configurations are set.\n");

    int total_size = 81;
    vvm_index vvms_index[total_size];
    int count = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i <= NX-2; i++) {
            for (int j = 2; j <= NY-2; j++) {
                if (p == 1 && (NX/2-30 <= i && i <= NX/2+30) && j == NY/2) {
                    vvms_index[count] = {p, i, j};
                    count++;
                }
                if (p == 1 && i == NX/2 && (NY/2-10 <= j && j <= NY/2+10 && j != NY/2)) {
                    vvms_index[count] = {p, i, j};
                    count++;
                }

                // if (p == 1 && (i >= NX/2-10 && i <= NX/2+10) && (j >= NX/2-10 && j <= NX/2+10)) {
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
    
    vvm**** vvms = allocate_and_initialize();
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
    double temp_csswm = TIMEEND / DT, temp_vvm = TIMEEND / config_vvms[1][NX/2][NY/2]->dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    CSSWM::Outputs::create_all_directory();
    // create Q_all directory
    CSSWM::Outputs::create_directory(std::string(OUTPUTPATH) + "Q_all/");

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
    
    double coupling_csswm_param = 2.;
    double coupling_vvm_param = 4.55;

    double thm_mean = 0.;
    double th_mean = 0.;
    double Q_all[6][NX][NY];
    #if defined(AB2_COUPLE)
        double q_all[2][6][NX][NY];
    #else
        double q_all[6][NX][NY];
    #endif
    // initialize Q_all, q_all
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                Q_all[p][i][j] = 0.;
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
        double time_csswm = n_csswm * DT;
        printf("n_vvm: %d, time_vvm: %f, n_csswm: %d,  time_csswm: %f\n", vvms[1][47][47]->step, time_vvm, n_csswm, time_csswm);

        if (time_vvm == time_csswm) {
            if (n_csswm % OUTPUTINTERVAL == 0 || n_csswm == TIMEEND-1 || n_csswm == TIMEEND-2) {
                #ifdef NCOUTPUT
                    CSSWM::Outputs::huv_nc(n_csswm, model_csswm);
                #endif
            }

            n_csswm++;

            CSSWM::Iteration::ph_pt_4(model_csswm);
            #ifndef Advection
                CSSWM::Iteration::pu_pt_4(model_csswm);
                CSSWM::Iteration::pv_pt_4(model_csswm);
            #endif

            // Exchange information here
            #ifdef _OPENMP
            #pragma omp parallel for reduction(+:thm_mean)
            #endif
            for (int size = 0; size < total_size; size++) {
                int p = vvms_index[size].p;
                int i = vvms_index[size].i;
                int j = vvms_index[size].j;
                thm_mean = 0.;
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        thm_mean += vvms[p][i][j]->thm[i_vvm][k_vvm];
                    }
                }
                thm_mean /= ((vvm_nx-2) * (vvm_nz-2));

                Q_all[p][i][j] = (exchange_coeff * thm_mean - model_csswm.csswm[p].hm[i][j]) / D2T;

                model_csswm.csswm[p].hp[i][j] += coupling_csswm_param * Q_all[p][i][j] * D2T;
            }
            #ifdef _OPENMP
            #pragma omp barrier
            #endif

            output_Qall(OUTPUTPATH + (std::string) "Q_all/", n_csswm, Q_all);

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

            
            #if defined(TIMEFILTER)
                CSSWM::NumericalProcess::timeFilterAll(model_csswm);
            #endif

            // next step
            CSSWM::Iteration::nextTimeStep(model_csswm);
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

            #if defined(AB2_COUPLE)
                output_qall(path + "vvm/q_all/", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, q_all[1]);
            #else
                output_qall(path + "vvm/q_all/", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, q_all);
            #endif
        }

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;
            vvms[p][i][j]->step++;
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

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
            #if defined(AB2_COUPLE)
                q_all[(vvms[p][i][j]->step+1)%2][p][i][j] = coupling_vvm_param * (model_csswm.csswm[p].hp[i][j] / exchange_coeff - th_mean) / DT;
                if (vvms[p][i][j]->step == 1) q_all[1][p][i][j] = q_all[0][p][i][j];
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        vvms[p][i][j]->thp[i_vvm][k_vvm] += 1.5*vvms[p][i][j]->dt*q_all[(vvms[p][i][j]->step+1)%2][p][i][j] - 0.5*vvms[p][i][j]->dt*q_all[vvms[p][i][j]->step%2][p][i][j];
                    }
                }
            #else
                q_all[p][i][j] = coupling_vvm_param * (model_csswm.csswm[p].hp[i][j] / exchange_coeff - th_mean) / D2T;
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        vvms[p][i][j]->thp[i_vvm][k_vvm] += vvms[p][i][j]->dt * q_all[p][i][j];
                    }
                }
            #endif
            
            // output q all nc
            vvm::Iteration::nextTimeStep(*vvms[p][i][j]);
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
    }

    deallocate_config(config_vvms);
    deallocate(vvms);
    return 0;
}

void output_qall(std::string dir, int n, double q[6][NX][NY]) {
    NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);
    // Create netCDF dimensions
    NcDim p = dataFile.addDim("p", 6);
    NcDim xDim = dataFile.addDim("x", NX);
    NcDim yDim = dataFile.addDim("y", NY);
    NcDim lonDim = dataFile.addDim("lon", NX);
    NcDim latDim = dataFile.addDim("lat", NY);

    std::vector<NcDim> xyDim, lonlatDim;
    xyDim.push_back(p);
    xyDim.push_back(xDim);
    xyDim.push_back(yDim);

    NcVar q_all = dataFile.addVar("q", ncDouble, xyDim);

    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(NX);
    countp.push_back(NY);

    for (int p = 0; p < 6; p++) {
        startp[0] = p;
        q_all.putVar(startp, countp, q[p]);
    }
    return;
}

void output_Qall(std::string dir,int n, double Q[6][NX][NY]) {
    NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);       
    // Create netCDF dimensions
    NcDim p = dataFile.addDim("p", 6);
    NcDim xDim = dataFile.addDim("x", NX);
    NcDim yDim = dataFile.addDim("y", NY);
    NcDim lonDim = dataFile.addDim("lon", NX);
    NcDim latDim = dataFile.addDim("lat", NY);

    std::vector<NcDim> xyDim, lonlatDim;
    xyDim.push_back(p);
    xyDim.push_back(xDim);
    xyDim.push_back(yDim);

    NcVar q_all = dataFile.addVar("q", ncDouble, xyDim);

    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(NX);
    countp.push_back(NY);

    for (int p = 0; p < 6; p++) {
        startp[0] = p;
        q_all.putVar(startp, countp, Q[p]);
    }
    return;
}
