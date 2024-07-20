#include "../2DVVM/src/Declare.hpp"
#include "../CSSWM/src/construction.hpp"
#include "../CSSWM/src/define.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <netcdf>
#include <mpi.h>

#if defined(PETSC)
    #include <petsc.h>
#endif

using namespace netCDF;


// CASE0: Nothing, CASE1:Bubble
Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double Kx, double Kz) {
    return Config_VVM(3.0, 200.0, 200.0, 100000, 20000, 1500000.0, 
                      10000, path, 10, 
                      Kx, Kz, 0.01, 0.0, 0.0, 1E-22, 9.80665, 1003.5, 716.5, 287.0, 
                      2.5E6, 1E5, 96500.0, addforcingtime, CASE);
}

// void output_qall(std::string dir,int n, double q[6][nx][ny]);
// void output_Qall(std::string dir,int n, double Q[6][nx][ny]);

#define MASTER_RANK 0

typedef struct {
    int p, i, j;
} vvm_index;

typedef struct {
    int p, i, j;
    double value;
} data_send;

MPI_Datatype create_data_send_type() {
    MPI_Datatype data_send_type;
    int lengths[2] = {3, 1}; // 3 ints, 1 double
    const MPI_Aint displacements[2] = {offsetof(data_send, p), offsetof(data_send, value)};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};

    MPI_Type_create_struct(2, lengths, displacements, types, &data_send_type);
    MPI_Type_commit(&data_send_type);

    return data_send_type;
}

void validate_index(vvm_index index, int nx_internal, int ny_internal) {
    if (index.p < 0 || index.p >= 6) {
        printf("Error: p=%d is out of range.\n", index.p);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (index.i < 2 || index.i >= 2 + nx_internal) {
        printf("Error: i=%d is out of range.\n", index.i);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (index.j < 2 || index.j >= 2 + ny_internal) {
        printf("Error: j=%d is out of range.\n", index.j);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int is_in_index_array(vvm_index value, vvm_index *array, int size) {
    for (int i = 0; i < size; i++) {
        if (array[i].p == value.p && array[i].i == value.i && array[i].j == value.j) {
            return 1;
        }
    }
    return 0;
}

int compute_global_index(vvm_index index, int nx_internal, int ny_internal) {
    return index.p * nx_internal * ny_internal + (index.i - 2) * ny_internal + (index.j - 2);
}

int main(int argc, char **argv) {
    CSSWM *model_csswm = nullptr;

    double exchange_coeff = 287. / 9.80665;
    double coupling_csswm_param = 2.;
    double coupling_vvm_param = 4.7;

    int total_size = 81;
    vvm_index valid_indices[total_size];
    int count = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i <= NX-3; i++) {
            for (int j = 2; j <= NY-3; j++) {
                // if (p == 1 && j == NY/2) {
                if (p == 1 && (NX/2-30 <= i && i <= NX/2+30) && j == NY/2) {
                    valid_indices[count] = {p, i, j};
                    count++;
                }

                if (p == 1 && i == NX/2 && (NY/2-10 <= j && j <= NY/2+10 && j != NY/2)) {
                    valid_indices[count] = {p, i, j};
                    count++;
                }
            }
        }
    }
    int valid_indices_size = sizeof(valid_indices) / sizeof(valid_indices[0]);
    if (valid_indices_size != count) {
        printf("Error: Incorrect valid_indices_size. Expected %d, got %d.\n", total_size, count);
        return 1;
    }


    MPI_Init(&argc, &argv);
    std::string path = "/data/Aaron/TMIF/newTur2_dt_4.7_cloud_2_csswm_1E6diff_60p1/";
    std::string path_vvm;

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype data_send_type = create_data_send_type();

    int tasks_per_core = valid_indices_size / size;
    int extra_tasks = valid_indices_size % size;

    int start, end;
    if (rank < extra_tasks) {
        // Some ranks handle one extra task
        start = rank * (tasks_per_core + 1);
        end = start + tasks_per_core;
    }
    else {
        // Other ranks handle standard tasks
        start = rank * tasks_per_core + extra_tasks;
        end = start + tasks_per_core - 1;
    }

    printf("rank=%d, start=%d, end=%d, load=%d\n", rank, start, end, end-start+1);

    // Gather info for master
    int* gather_counts = new int[size];
    int* displs_gather = new int[size];

    // Master process: prepare to receive the data
    for (int i = 0; i < size; ++i) {
        if (i < extra_tasks) gather_counts[i] = tasks_per_core + 1;
        else gather_counts[i] = tasks_per_core;
    }
    // Calculate displacements
    int total_gather = 0;
    for (int i = 0; i < size; ++i) {
        displs_gather[i] = total_gather;
        total_gather += gather_counts[i];
    }
    if (rank == MASTER_RANK) {
        if (gather_counts[0] != end-start+1) {
            printf("Error: Incorrect gather_counts[0]. Expected %d, got %d.\n", end-start+1, gather_counts[0]);
            return 1;
        }
    }
    
    // Scatter info for all
    int* scatter_counts = new int[size];
    int* displs_scatter = new int[size];

    for (int i = 0; i < size; i++) {
        if (i < extra_tasks) scatter_counts[i] = tasks_per_core + 1;
        else scatter_counts[i] = tasks_per_core;
    }

    int total_scatter = 0;
    for (int i = 0; i < size; ++i) {
        displs_scatter[i] = total_scatter;
        total_scatter += scatter_counts[i];
    }
    data_send scattered_data[total_scatter];
    data_send received_data[scatter_counts[rank]];
    double q_AB2[scatter_counts[rank]][2];

    
    if (rank == MASTER_RANK) {
        model_csswm = new CSSWM();
        CSSWM::Init::Init2d(*model_csswm);

        CSSWM::Outputs::create_all_directory();
        CSSWM::Outputs::create_directory(path + "csswm/Q_all/");
        #ifdef NCOUTPUT
            CSSWM::Outputs::grid_nc(*model_csswm);
        #endif
        #ifdef TXTOUTPUT
            CSSWM::Outputs::grid(model_csswm);
        #endif

        vvm::Output::create_directory(path + "vvm/q_all/");
    }

    vvm **vvms = new vvm*[end-start+1];
    for (int v = start; v <= end; v++) {
        if (v >= valid_indices_size) break;
        int local_v = v - start;

        vvm_index current_index = valid_indices[v];
        printf("Rank %d: Processing p=%d, i=%d, j=%d\n", rank, current_index.p, current_index.i, current_index.j);
        path_vvm = path + "vvm/" + std::to_string(current_index.p) + "_" + std::to_string(current_index.i) + "_" + std::to_string(current_index.j) + "/";

        if (current_index.p == 1 && current_index.i == NX/2 && (NY/2-10 <= current_index.j && current_index.j <= NY/2+10)) {
            vvms[local_v] = new vvm(createConfig(path_vvm, 10, 1, 200, 200));
        }
        else {
            vvms[local_v] = new vvm(createConfig(path_vvm, 10, 0, 70, 70));
        }
        vvm::Output::create_all_directory(*vvms[local_v]);


        vvm::Init::Init1d(*vvms[local_v]);
        vvm::Init::Init2d(*vvms[local_v]);
        #ifndef PETSC
            vvm::PoissonSolver::InitPoissonMatrix(*vvms[local_v]);
        #endif
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int n_csswm = 0;
    double temp_csswm = TIMEEND / DT, temp_vvm = TIMEEND / vvms[0]->dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    while (vvms[0]->step < nmax_vvm || n_csswm < nmax_csswm) {
        double time_vvm = vvms[0]->step * vvms[0]->dt;
        double time_csswm = n_csswm * DT;
        if (rank == MASTER_RANK) printf("n_vvm: %d, time_vvm: %f, n_csswm: %d,  time_csswm: %f\n", vvms[0]->step, time_vvm, n_csswm, time_csswm);


        if (time_vvm == time_csswm) {
            MPI_Barrier(MPI_COMM_WORLD);

            if (rank == MASTER_RANK) {
                if (n_csswm % OUTPUTINTERVAL == 0 || n_csswm == TIMEEND-1 || n_csswm == TIMEEND-2) {
                    #ifdef NCOUTPUT
                        CSSWM::Outputs::huv_nc(n_csswm, *model_csswm);
                    #endif
                }

                n_csswm++;

                CSSWM::Iteration::ph_pt_4(*model_csswm);
                #ifndef Advection
                    CSSWM::Iteration::pu_pt_4(*model_csswm);
                    CSSWM::Iteration::pv_pt_4(*model_csswm);
                #endif

                // TODO: Exchange information here
                data_send* gathered_data = new data_send[total_gather];
                MPI_Request* requests = new MPI_Request[size - 1];
                for (int i = 1; i < size; ++i) {
                    MPI_Irecv(&gathered_data[displs_gather[i]], gather_counts[i], data_send_type, i, 0, MPI_COMM_WORLD, &requests[i - 1]);
                }

                for (int v = start; v <= end; v++) {
                    if (v >= valid_indices_size) break;
                    vvm_index current_index = valid_indices[v];
                    int local_v = v - start;

                    double thm_mean = 0.;
                    for (int k_vvm = 1; k_vvm <= vvms[local_v]->nz-2; k_vvm++) {
                        for (int i_vvm = 1; i_vvm <= vvms[local_v]->nx-2; i_vvm++) {
                            thm_mean += vvms[local_v]->thm[i_vvm][k_vvm];
                        }
                    }
                    thm_mean /= ((vvms[local_v]->nx-2) * (vvms[local_v]->nz-2));

                    gathered_data[local_v] = {current_index.p, current_index.i, current_index.j, thm_mean};
                }

                for (int i = 1; i < size; ++i) {
                    MPI_Wait(&requests[i - 1], MPI_STATUS_IGNORE);
                }
                
                // MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                //             gathered_data, gather_counts, displs_gather, data_send_type, MASTER_RANK, MPI_COMM_WORLD);


                // for (int i = 0; i < total_gather; ++i) {
                //     std::cout << "Received value " << gathered_data[i].value
                //             << " with indices (" << gathered_data[i].p << ", "
                //             << gathered_data[i].i << ", " << gathered_data[i].j << ")" << std::endl;
                // }

                // Exchange here
                double Q = 0.;
                for (int i = 0; i < total_gather; i++) {
                    Q = (exchange_coeff * gathered_data[i].value - model_csswm->csswm[gathered_data[i].p].hm[gathered_data[i].i][gathered_data[i].j]) / D2T;
                    model_csswm->csswm[gathered_data[i].p].hp[gathered_data[i].i][gathered_data[i].j] += coupling_csswm_param * Q * D2T;
                }


                delete[] gathered_data;
                delete[] requests;

                // Boundary exchange and interpolation
                model_csswm->BP_h(*model_csswm);
                #ifndef Advection
                    model_csswm->BP_wind_interpolation2(*model_csswm);
                #endif

                #if defined(DIFFUSION)
                    CSSWM::NumericalProcess::DiffusionAll(*model_csswm);
                #endif

                model_csswm->BP_h(*model_csswm);
                model_csswm->BP_wind_interpolation2(*model_csswm);
                
                #if defined(TIMEFILTER)
                    CSSWM::NumericalProcess::timeFilterAll(*model_csswm);
                #endif

                // next step
                CSSWM::Iteration::nextTimeStep(*model_csswm);
            }
            else {
                n_csswm++;

                // Worker processes: send data
                int send_size = end - start + 1;
                data_send* send_data = new data_send[send_size];

                for (int v = start; v <= end; v++) {
                    if (v >= valid_indices_size) break;
                    vvm_index current_index = valid_indices[v];

                    int local_v = v - start;
                    double thm_mean = 0.;
                    for (int k_vvm = 1; k_vvm <= vvms[local_v]->nz-2; k_vvm++) {
                        for (int i_vvm = 1; i_vvm <= vvms[local_v]->nx-2; i_vvm++) {
                            thm_mean += vvms[local_v]->thm[i_vvm][k_vvm];
                        }
                    }
                    thm_mean /= ((vvms[local_v]->nx-2) * (vvms[local_v]->nz-2));

                    send_data[local_v] = {current_index.p, current_index.i, current_index.j, thm_mean};
                }

                MPI_Request send_request;
                MPI_Isend(send_data, send_size, data_send_type, MASTER_RANK, 0, MPI_COMM_WORLD, &send_request);
                MPI_Wait(&send_request, MPI_STATUS_IGNORE);

                // MPI_Barrier(MPI_COMM_WORLD);

                // MPI_Gatherv(send_data, send_size, data_send_type,
                //             nullptr, nullptr, nullptr, data_send_type, MASTER_RANK, MPI_COMM_WORLD);


                // Clean up allocated memory
                delete[] send_data;
            }
        }

        if (rank == MASTER_RANK) {
            for (int i = 0; i < total_scatter; i++) {
                scattered_data[i] = {valid_indices[i].p, valid_indices[i].i, valid_indices[i].j, model_csswm->csswm[valid_indices[i].p].hp[valid_indices[i].i][valid_indices[i].j]};
            }

            for (int i = 1; i < size; i++) {
                MPI_Send(&scattered_data[displs_scatter[i]], scatter_counts[i], data_send_type, i, 0, MPI_COMM_WORLD);
            }

            for (int i = 0; i < scatter_counts[0]; i++) {
                received_data[i] = {valid_indices[i].p, valid_indices[i].i, valid_indices[i].j, model_csswm->csswm[valid_indices[i].p].hp[valid_indices[i].i][valid_indices[i].j]};
            }
        }
        else {
            MPI_Recv(received_data, scatter_counts[rank], data_send_type, MASTER_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("Rank %d received: ", rank);
            // for (int i = 0; i < scatter_counts[rank]; ++i) {
            //     printf("(%d, %d, %d, %f) ", received_data[i].p, received_data[i].i, received_data[i].j, received_data[i].value);
            // }
            // printf("\n");
        }



        for (int v = start; v <= end; v++) {
            if (v >= valid_indices_size) break;
            int local_v = v - start;
            vvm_index current_index = valid_indices[v];

            if (vvms[local_v]->step % vvms[local_v]->OUTPUTSTEP == 0) {
                #if defined(OUTPUTTXT)
                    vvm::Output::outputalltxt(vvms[local_v]->step, *vvms[local_v]);
                #endif

                #if defined(OUTPUTNC)
                    vvm::Output::output_nc(vvms[local_v]->step, *vvms[local_v]);
                #endif
            }

            // TODO: Output Exchange information here
            // #if defined(AB2)
            //     output_qall(path + "vvm/q_all/", vvms[local_v]->step, q_all[1]);
            // #else
            //     output_qall(path + "vvm/q_all/", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, q_all);
            // #endif

            vvms[local_v]->step++;

            vvm::Iteration::pzeta_pt(*vvms[local_v]);
            vvm::Iteration::pth_pt(*vvms[local_v]);
            #if defined(WATER)
                vvm::Iteration::pqv_pt(*vvms[local_v]);
                vvm::Iteration::pqc_pt(*vvms[local_v]);
                vvm::Iteration::pqr_pt(*vvms[local_v]);

                if (vvms[local_v]->step * vvms[local_v]->dt <= vvms[local_v]->addforcingtime) vvms[local_v]->status_for_adding_forcing = true;
                else vvms[local_v]->status_for_adding_forcing = false;

                // Generate new random th perturbation for tropical forcing case
                if (vvms[local_v]->status_for_adding_forcing == true) {
                    vvm::Init::RandomPerturbation(*vvms[local_v], vvms[local_v]->step);
                }
                vvm::AddForcing(*vvms[local_v]);
            #endif
            vvm::BoundaryProcess2D_all(*vvms[local_v]);

            vvm::PoissonSolver::pubarTop_pt(*vvms[local_v]);
            vvm::PoissonSolver::cal_w(*vvms[local_v], current_index.p, current_index.i, current_index.j);
            vvm::PoissonSolver::cal_u(*vvms[local_v]);
            
            vvm::Iteration::updateMean(*vvms[local_v]);
            vvm::Turbulence::RKM_RKH(*vvms[local_v]);

            #if defined(WATER)
                vvm::MicroPhysics::autoconversion(*vvms[local_v]);
                vvm::MicroPhysics::accretion(*vvms[local_v]);
                vvm::MicroPhysics::evaporation(*vvms[local_v]);
                vvm::MicroPhysics::condensation(*vvms[local_v]); // saturation adjustment

                // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
                vvm::MicroPhysics::NegativeValueProcess(vvms[local_v]->qvp, vvms[local_v]->nx, vvms[local_v]->nz);
                vvm::MicroPhysics::NegativeValueProcess(vvms[local_v]->qcp, vvms[local_v]->nx, vvms[local_v]->nz);
                vvm::MicroPhysics::NegativeValueProcess(vvms[local_v]->qrp, vvms[local_v]->nx, vvms[local_v]->nz);
            #endif

            vvm::BoundaryProcess2D_all(*vvms[local_v]);

            #if defined(TIMEFILTER) && !defined(AB2)
                vvm::NumericalProcess::timeFilterAll(*vvms[local_v]);
            #endif

        }
        MPI_Barrier(MPI_COMM_WORLD);

        // TODO: Exchange information here
        for (int v = start; v <= end; v++) {
            if (v >= valid_indices_size) break;
            int local_v = v - start;

            // Exchange here
            double th_mean = 0.;
            for (int k_vvm = 1; k_vvm <= vvms[local_v]->nz-2; k_vvm++) {
                for (int i_vvm = 1; i_vvm <= vvms[local_v]->nx-2; i_vvm++) {
                    th_mean += vvms[local_v]->th[i_vvm][k_vvm];
                }
            }
            th_mean /= ((vvms[local_v]->nx-2) * (vvms[local_v]->nz-2));

            q_AB2[local_v][(vvms[local_v]->step+1)%2] = coupling_vvm_param * (received_data[local_v].value / exchange_coeff - th_mean) / DT;
            if (vvms[local_v]->step == 1) q_AB2[local_v][1] = q_AB2[local_v][0];
            for (int k_vvm = 1; k_vvm <= vvms[local_v]->nz-2; k_vvm++) {
                for (int i_vvm = 1; i_vvm <= vvms[local_v]->nx-2; i_vvm++) {
                    vvms[local_v]->thp[i_vvm][k_vvm] += 1.5*vvms[local_v]->dt*q_AB2[local_v][(vvms[local_v]->step+1)%2] - 0.5*vvms[local_v]->dt*q_AB2[local_v][vvms[local_v]->step%2];
                }
            }

            vvm::Iteration::nextTimeStep(*vvms[local_v]);
        }
    }


    for (int v = 0; v < end-start+1; v++) {
        if (vvms[v] != nullptr) {
            delete vvms[v];
        }
    }
    delete[] vvms;
    delete[] gather_counts;
    delete[] displs_gather;

    MPI_Type_free(&data_send_type);
    MPI_Finalize();

    return 0;
}