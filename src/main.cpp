#include "allocate_csswm_vvms.hpp"
#include "reading_config.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif
#if defined(NCOUTPUT)
    #include <netcdf>
#endif
#include <cstdio>
#include <mpi.h>

using namespace netCDF;

#define MASTER_RANK 0



typedef struct {
    int p, i, j;
    double value;
} data_send;

MPI_Datatype create_data_send_type() {
    MPI_Datatype data_send_type;
    int lengths[2] = {3, 1};
    const MPI_Aint displacements[2] = {offsetof(data_send, p), offsetof(data_send, value)};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};

    MPI_Type_create_struct(2, lengths, displacements, types, &data_send_type);
    MPI_Type_commit(&data_send_type);
    return data_send_type;
}

#if defined(P3_MICROPHY)
extern "C" {
    void __microphy_p3_MOD_p3_init(
        char* lookup_file_dir, int* nCat, bool* trplMomI, bool* liqfrac,
        char* model, int* stat, bool* abort_on_err, bool* dowr, size_t lookup_file_dir_len, size_t model_name_len
    );
}


extern "C" {
    void __microphy_p3_MOD_p3_main(
        double* qc, double* nc, double* qr, double* nr, 
        double* th_old, double* th, double* qv_old, double* qv, 
        double* dt, double* qitot, double* qirim, double* qiliq, 
        double* ni, double* birim, double* zi, double* ssat, 
        double* w, double* p, double* dz, int* itt, 
        double* precip_liq, double* precip_sol, 
        int* one, int* ncols, int* one2, int* nz, int* nCat, 
        double* diag_ze, double* diag_effc, double* diag_effi, 
        double* diag_vmi, double* diag_di, double* diag_rhoi, 
        int* n_diag_2d, double* diag_2d, int* n_diag_3d, double* diag_3d, 
        bool* log_predictNc, char* model, 
        double* clbfact_dep, double* clbfact_sub, 
        bool* debug_on, bool* scpf_on, double* scpf_pfrac, 
        double* scpf_resfact, double* cldfrac, 
        bool* trplMomI, bool* liqfrac,
        double* , double*, double*, double*, double*,
        double* , double*, double*, double*, double*,
        double* , double*, double*, double*, double*, size_t model_name_len
    );
}
#endif

void output_forcing(std::string dir, int n, double q[6][NX][NY]);
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype data_send_type = create_data_send_type();

    // This heating weight follows the Q1 heating profile for the data in 2DVVM/input/init.txt
    // Couple 12 km
    double heating_weight[62] = {
        0.,
        0.0012225044493978194, 0.004229087830248994, 0.00721407219255957, 0.009589964232025298, 0.01192265823440983, 0.01382337186598241, 0.015335303163824239, 0.016847234461666062, 0.017840789314533548, 0.018661552019076256, 
        0.019439116686537767, 0.020130285279836888, 0.02060546368773003, 0.021080642095623176, 0.021555820503516322, 0.02198780087432827, 0.022290187133896636, 0.0226357714305462, 0.022981355727195757, 0.023283741986764125, 
        0.023456534135088903, 0.023586128246332487, 0.023715722357576074, 0.02384531646881966, 0.02393171254298205, 0.023758920394657268, 0.023629326283413684, 0.0234997321721701, 0.023370138060926512, 0.02324054394968293, 
        0.022938157690114563, 0.0226357714305462, 0.022376583208059027, 0.022117394985571855, 0.021815008726003494, 0.02146942442935393, 0.02103744405854198, 0.02060546368773003, 0.02017348331691808, 0.019784700983187326, 
        0.019352720612375373, 0.018575155944913865, 0.017797591277452354, 0.017063224647072037, 0.01628565997961053, 0.015551293349230213, 0.014730530644687507, 0.013736975791820021, 0.012743420938952534, 0.01174986608608505, 
        0.010681696441895499, 0.009613526797705949, 0.0085453571535164, 0.00747718750932685, 0.0064090178651373, 0.005340848220947749, 0.0042726785767582005, 0.003204508932568649, 0.0021363392883791002, 0.0010681696441895501,
        0.
    };

    Eigen::setNbThreads(1);
    
    // Declare CSSWM
    CSSWM *model_csswm = nullptr;
    if (rank == MASTER_RANK) model_csswm = new CSSWM();

    // Read configuration file
    std::map<std::string, std::string> config = read_config("../config.txt");
    if (rank == MASTER_RANK) {
        std::cout << "Configuration key-value pairs:" << std::endl;
        for (const auto& pair : config) {
            std::cout << pair.first << " = " << pair.second << std::endl;
        }
    }

    std::string path = config["OUTPUTPATH"] + "/";
    int seed = std::stoi(config["SEED"]);
    double Couple_time = std::stod(config["COUPLETIME"]);
    std::vector<vvm_index> Bubbles_p_i_j = parse_int_tuples(config["Bubble_p_i_j"]);
    std::vector<vvm_index> NotBubbles_p_i_j = parse_int_tuples(config["NotBubble_p_i_j"]);
    int Bubble_case = std::stoi(config["BubbleCase"]);
    double csswm_gravity = std::stod(config["CSSWM_GRAVITY"]);

    double csswm_dt = std::stod(config["CSSWM_DT"]);
    double csswm_timeend = std::stod(config["CSSWM_TIMEEND"]);
    int csswm_outputstep = std::stoi(config["CSSWM_OUTPUTSTEP"]);
    double csswm_diffusion_kx = std::stod(config["CSSWM_DIFFUSION_KX"]);
    double csswm_diffusion_ky = std::stod(config["CSSWM_DIFFUSION_KY"]);
    double csswm_diffusion_ts = std::stod(config["CSSWM_DIFFUSION_TS"]);
    double csswm_addforcingtime = std::stod(config["CSSWM_ADDFORCING_TIME"]);
    double csswm_h_nudge_time = std::stod(config["CSSWM_H_NUDGE_TIME"]);

    double vvm_xrange = std::stod(config["VVM_XRANGE"]);
    double vvm_zrange = std::stod(config["VVM_ZRANGE"]);
    double vvm_dx = std::stod(config["VVM_DX"]);
    double vvm_dz = std::stod(config["VVM_DZ"]);
    double vvm_dt = std::stod(config["VVM_DT"]);
    double vvm_timeend = std::stod(config["VVM_TIMEEND"]);
    int vvm_outputstep = std::stoi(config["VVM_OUTPUTSTEP"]);
    double vvm_moisture_nudge_time = std::stod(config["VVM_MOISTURE_NUDGE_TIME"]);


    if (rank == MASTER_RANK) {
        model_csswm->output_path = path + "csswm/";
        model_csswm->gravity = csswm_gravity;
        model_csswm->dt = csswm_dt;
        model_csswm->timeend = csswm_timeend;
        model_csswm->outputstep = csswm_outputstep;
        model_csswm->diffusion_kx = csswm_diffusion_kx;
        model_csswm->diffusion_ky = csswm_diffusion_ky;
        model_csswm->diffusion_ts = csswm_diffusion_ts;
        model_csswm->addforcingtime = csswm_addforcingtime;
        model_csswm->csswm_h_nudge_time = csswm_h_nudge_time;
        
    }

    std::vector<vvm_index> vvms_index;
    for (const auto& bubble : Bubbles_p_i_j) {
        vvms_index.push_back(bubble);
    }
    for (const auto& notbubble : NotBubbles_p_i_j) {
        vvms_index.push_back(notbubble);
    }
    int total_vvm_size = vvms_index.size();

    // Distribute work
    int tasks_per_core = total_vvm_size / size;
    int extra_tasks = total_vvm_size % size;
    int start = rank * tasks_per_core + (rank < extra_tasks ? rank : extra_tasks);
    int end = start + tasks_per_core + (rank < extra_tasks ? 1 : 0) - 1;
    if (end >= total_vvm_size) end = total_vvm_size - 1;

    // Gather/scatter setup
    int* gather_counts = new int[size];
    int* displs_gather = new int[size];
    int* scatter_counts = new int[size];
    int* displs_scatter = new int[size];
    
    for (int i = 0; i < size; i++) {
        gather_counts[i] = scatter_counts[i] = (i < extra_tasks) ? tasks_per_core + 1 : tasks_per_core;
        displs_gather[i] = (i == 0) ? 0 : displs_gather[i-1] + gather_counts[i-1];
        displs_scatter[i] = (i == 0) ? 0 : displs_scatter[i-1] + scatter_counts[i-1];
    }

    data_send* scattered_data = new data_send[total_vvm_size];
    data_send* received_data = new data_send[scatter_counts[rank]];

    // Master initializes CSSWM
    if (rank == MASTER_RANK) {
        CSSWM::Init::Init2d(*model_csswm);
        CSSWM::Outputs::create_all_directory(*model_csswm);
        CSSWM::Outputs::create_directory(model_csswm->output_path + "Q_all/");
        CSSWM::Outputs::create_directory(model_csswm->output_path + "q_all/");
        #ifdef NCOUTPUT
            CSSWM::Outputs::grid_nc(*model_csswm);
        #endif
        #ifdef TXTOUTPUT
            CSSWM::Outputs::grid(model_csswm);
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Local vvm declaration
    int local_size = end - start + 1;
    Config_VVM **config_vvms = new Config_VVM*[local_size];
    vvm** local_vvms = new vvm*[local_size];

    // Local VVM initialization
    std::string path_vvm;
    int vvm_nx = 0, vvm_nz = 0;
    
    for (int v = start; v <= end; v++) {
        vvm_index current_index = vvms_index[v];
        int p = vvms_index[v].p;
        int i = vvms_index[v].i;
        int j = vvms_index[v].j;
        int local_v = v - start;

        path_vvm = path + "vvm/" + std::to_string(p) + "_" + std::to_string(i) + "_" + std::to_string(j) + "/";
        bool is_bubble = std::find(Bubbles_p_i_j.begin(), Bubbles_p_i_j.end(), vvms_index[v]) != Bubbles_p_i_j.end();
        
        printf("Rank %d: Processing p=%d, i=%d, j=%d. Is bubble: %d\n", rank, p, i, j, is_bubble);
        config_vvms[local_v] = new Config_VVM(createConfig(path_vvm, 10, 
                                                is_bubble ? Bubble_case : 0, vvm_xrange, vvm_zrange, vvm_dx, vvm_dz, 
                                                vvm_dt, vvm_timeend, vvm_outputstep, vvm_moisture_nudge_time));

        local_vvms[local_v] = new vvm(*config_vvms[local_v]);


        #if defined(LOADFROMPREVIOUSFILE)
            vvm::Init::LoadFromPreviousFile(*vvms[p][i][j]);
        #elif defined(LOAD2DINIT)
            vvm::Init::Load2DInit(*vvms[p][i][j]);
        #else
            vvm::Init::Init1d(*local_vvms[local_v]);
            vvm::Init::Init2d(*local_vvms[local_v]);
            #ifndef PETSC
                vvm::PoissonSolver::InitPoissonMatrix(*local_vvms[local_v]);
            #endif
        #endif
        vvm::Output::create_all_directory(*local_vvms[local_v]);
        vvm::Output::create_directory(path + "vvm/q_all/");

        vvm_nx = local_vvms[local_v]->nx;
        vvm_nz = local_vvms[local_v]->nz;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank: %d. Configurations are set and VVMs are declared and initialized.\n", rank);



    // TODO: Flexible ctl output
    // Copy grads ctl file to the output directory
    if (rank == MASTER_RANK) {
        std::string src1 = "../CSSWM/scripts/csswm.ctl";
        std::string des1 = model_csswm->output_path + "nc/csswm.ctl";

        // Construct the command
        std::string cmd1 = "cp " + src1 + " " + des1;
        
        // Execute the command
        system(cmd1.c_str());
    }

    std::string src2 = "../2DVVM/scripts/vvm.ctl";
    std::string des2;
    std::string cmd2;

    for (int v = start; v <= end; v++) {
        int local_v = v - start;
        des2 = local_vvms[local_v]->outputpath + "nc/vvm.ctl";        
        cmd2 = "cp " + src2 + " " + des2;
        system(cmd2.c_str());
    }

    double th_mean = 0.;
    double exchange_coeff = 0.;
    double coupling_vvm_param = 1.;
    int k_couple = -9999;
    k_couple = 12000. / local_vvms[0]->dz;
    double (*th_mean_all)[NX][NY] = nullptr;
    double (*Q_all)[NX][NY] = nullptr;
    double (*q_all)[NX][NY] = nullptr;
    if (rank == MASTER_RANK) {
        th_mean_all = new double[6][NX][NY]();
        Q_all = new double[6][NX][NY]();
        q_all = new double[6][NX][NY]();
    }

    // Calculate initial forcing and tune the profile according to the exchange coefficient.
    data_send* initial_send_data = new data_send[local_size];
    for (int v = start; v <= end; v++) {
        int local_v = v - start;
        int p = vvms_index[v].p;
        int i = vvms_index[v].i;
        int j = vvms_index[v].j;

        double th_mean = 0.;
        for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                th_mean += local_vvms[local_v]->th[i_vvm][k_vvm];
            }
        }
        th_mean /= ((vvm_nx-2) * k_couple);
        initial_send_data[local_v] = {p, i, j, th_mean};
    }
    data_send* initial_gathered_data = (rank == MASTER_RANK) ? new data_send[total_vvm_size] : nullptr;
    MPI_Gatherv(initial_send_data, end-start+1, data_send_type,
                initial_gathered_data, gather_counts, displs_gather, data_send_type,
                MASTER_RANK, MPI_COMM_WORLD);


    if (rank == MASTER_RANK) {
        for (int count = 0; count < total_vvm_size; count++) {
            int p = initial_gathered_data[count].p;
            int i = initial_gathered_data[count].i;
            int j = initial_gathered_data[count].j;
            th_mean_all[p][i][j] = initial_gathered_data[count].value;
        }

        // Initial adjustment using th_mean_all
        exchange_coeff = model_csswm->csswm[NotBubbles_p_i_j[0].p].h[NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j] / 
                        th_mean_all[NotBubbles_p_i_j[0].p][NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j];

        for (int v = 0; v < total_vvm_size; v++) {
            int p = vvms_index[v].p;
            int i = vvms_index[v].i;
            int j = vvms_index[v].j;
            double total_heating = (model_csswm->csswm[p].h[i][j] / exchange_coeff - th_mean_all[p][i][j]) * (vvm_nz-2);
            scattered_data[v] = {p, i, j, total_heating};
        }
        delete[] initial_gathered_data;
    }
    delete[] initial_send_data;


    MPI_Scatterv(scattered_data, scatter_counts, displs_scatter, data_send_type,
                 received_data, scatter_counts[rank], data_send_type,
                 MASTER_RANK, MPI_COMM_WORLD);
    
    // Broadcast the exchange coefficient to all ranks
    MPI_Bcast(&exchange_coeff, 1, MPI_DOUBLE, MASTER_RANK, MPI_COMM_WORLD);
    std::cout << "exchange_coeff" << exchange_coeff << std::endl;

    // Apply initial heating
    for (int v = start; v <= end; v++) {
        int local_v = v - start;

        double total_heating = received_data[local_v].value;
        for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
            double heating = total_heating * heating_weight[k_vvm];
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                local_vvms[local_v]->th[i_vvm][k_vvm] += heating;
                local_vvms[local_v]->thm[i_vvm][k_vvm] = local_vvms[local_v]->th[i_vvm][k_vvm];
            }
        }
    }
    // p3 initialization
    #if defined(P3_MICROPHY)

    vvm::P3::lookup_file_dir = strdup("../2DVVM/lookup_tables");
    __microphy_p3_MOD_p3_init(
        vvm::P3::lookup_file_dir, &vvm::P3::nCat, &vvm::P3::trplMomI, &vvm::P3::liqfrac,
        vvm::P3::model_name, &vvm::P3::stat, &vvm::P3::abort_on_err, &vvm::P3::dowr, strlen(vvm::P3::lookup_file_dir), strlen(vvm::P3::model_name)
    );

    for (int v = start; v <= end; v++) {
        int local_v = v - start;
        for (int k_vvm = 0; k_vvm <= vvm_nz-1; k_vvm++) {
            for (int i_vvm = 0; i_vvm <= vvm_nx-1; i_vvm++) {
                local_vvms[local_v]->dz_all[i_vvm][k_vvm] = local_vvms[local_v]->dz;
                local_vvms[local_v]->w_all[i_vvm][k_vvm] = 0.;
                local_vvms[local_v]->pb_all[i_vvm][k_vvm] = local_vvms[local_v]->pb[k_vvm];
                local_vvms[local_v]->zi_all[i_vvm][k_vvm] = 0.;
                local_vvms[local_v]->ssat_all[i_vvm][k_vvm] = 0.;
            }
        }
        for (int k_vvm = 0; k_vvm < vvm_nz; k_vvm++) {
            for (int i_vvm = 0; i_vvm < vvm_nx; i_vvm++) {
                local_vvms[local_v]->qiliqp[i_vvm][k_vvm] = 0.;
            }
        }
    }
    int one = 1;
    #endif
    MPI_Barrier(MPI_COMM_WORLD);

    double next_coupling_time = Couple_time;
    int n_csswm = 0;
    double temp_csswm = csswm_timeend / csswm_dt, temp_vvm = vvm_timeend / vvm_dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    while (true) {
        // Compute loop condition on all ranks for VVM and on master for CSSWM
        int vvm_continue = (local_vvms[0]->step < nmax_vvm);
        int csswm_continue = 1;
        if (rank == MASTER_RANK) csswm_continue = (model_csswm->step < nmax_csswm);
        MPI_Bcast(&csswm_continue, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
        if (!vvm_continue || !csswm_continue) break;

        double time_vvm = local_vvms[0]->step * local_vvms[0]->dt;
        double time_csswm = (rank == MASTER_RANK) ? model_csswm->step * DT : 0.0;

        if (rank == MASTER_RANK) {
            while (next_coupling_time != time_csswm) {
                printf("csswm_step: %d, csswm_time: %f\n", model_csswm->step, time_csswm);
                if (model_csswm->step % model_csswm->outputstep == 0 || 
                    model_csswm->step == model_csswm->timeend-1 || 
                    model_csswm->step == model_csswm->timeend-2) {
                    #ifdef NCOUTPUT
                        CSSWM::Outputs::huv_nc(model_csswm->step, *model_csswm);
                    #endif
                }

                #if defined(EquatorialWave)
                    if (model_csswm.step * model_csswm.dt >= model_csswm.addforcingtime) model_csswm.status_add_forcing = false;
                    else model_csswm.status_add_forcing = true;
                #endif

                // Prediction for CSSWM
                CSSWM::Iteration::ph_pt_4(*model_csswm);
                CSSWM::Iteration::pu_pt_4(*model_csswm);
                CSSWM::Iteration::pv_pt_4(*model_csswm);

                CSSWM::BP_h(*model_csswm);
                model_csswm->BP_wind_interpolation2(*model_csswm);

                #if defined(DIFFUSION)
                    CSSWM::NumericalProcess::DiffusionAll(*model_csswm);
                #endif

                model_csswm->BP_h(*model_csswm);
                model_csswm->BP_wind_interpolation2(*model_csswm);
                
                #if defined(TIMEFILTER) && !defined(AB2Time)
                    CSSWM::NumericalProcess::timeFilterAll(*model_csswm);
                #endif

                if (model_csswm->csswm_h_nudge_time != 0) {
                    CSSWM::NumericalProcess::NudgeH(*model_csswm);
                }

                CSSWM::Iteration::nextTimeStep(*model_csswm);
                model_csswm->step++;
                time_csswm = model_csswm->step * DT;
            }
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        MPI_Barrier(MPI_COMM_WORLD);

        // Master updates th_mean_all and computes forcing
        if (rank == MASTER_RANK) {
            for (int i = 0; i < total_vvm_size; i++) {
                int p = vvms_index[i].p;
                int x = vvms_index[i].i;
                int y = vvms_index[i].j;
                q_all[p][x][y] = (model_csswm->csswm[p].hp[x][y] / exchange_coeff - th_mean_all[p][x][y]) / Couple_time;
                // Q_all[p][x][y] = (exchange_coeff * th_mean_all[p][x][y] - model_csswm->csswm[p].h[x][y]) / Couple_time;
            }
            output_forcing(model_csswm->output_path + "q_all/", (int)next_coupling_time / Couple_time, q_all);
            // output_forcing(model_csswm->output_path + "Q_all/", (int)next_coupling_time / Couple_time, Q_all, rank);

            for (int i = 0; i < total_vvm_size; i++) {
                scattered_data[i] = {vvms_index[i].p, vvms_index[i].i, vvms_index[i].j, q_all[vvms_index[i].p][vvms_index[i].i][vvms_index[i].j]};
            }
        }

        MPI_Scatterv(scattered_data, scatter_counts, displs_scatter, data_send_type,
                     received_data, scatter_counts[rank], data_send_type,
                     MASTER_RANK, MPI_COMM_WORLD);


        while (time_vvm < next_coupling_time) {
            printf("VVM step: %d, time: %f\n", local_vvms[0]->step, time_vvm);
            
            for (int v = start; v <= end; v++) {
                int p = vvms_index[v].p;
                int i = vvms_index[v].i;
                int j = vvms_index[v].j;
                int local_v = v - start;

                if (local_vvms[local_v]->step % local_vvms[local_v]->OUTPUTSTEP == 0) {
                    #if defined(OUTPUTTXT)
                        vvm::Output::outputalltxt(vvms[p][i][j]->step, *vvms[p][i][j]);
                    #endif
                    #if defined(OUTPUTNC)
                        vvm::Output::output_nc(local_vvms[local_v]->step, *local_vvms[local_v]);
                    #endif
                }

                vvm::Iteration::pzeta_pt(*local_vvms[local_v]);
                vvm::Iteration::pth_pt(*local_vvms[local_v]);
                #if defined(WATER)
                    #if defined(KESSLER_MICROPHY)
                        vvm::Iteration::pqv_pt(*local_vvms[local_v]);
                        vvm::Iteration::pqc_pt(*local_vvms[local_v]);
                        vvm::Iteration::pqr_pt(*local_vvms[local_v]);
                    #endif

                    #if defined(P3_MICROPHY)
                        vvm::Iteration::pqmicrophy_pt(*local_vvms[local_v]);
                    #endif

                    if (local_vvms[local_v]->step * local_vvms[local_v]->dt <= local_vvms[local_v]->addforcingtime) local_vvms[local_v]->status_for_adding_forcing = true;
                    else local_vvms[local_v]->status_for_adding_forcing = false;


                    if (local_vvms[local_v]->status_for_adding_forcing) {
                        if (std::find(Bubbles_p_i_j.begin(), Bubbles_p_i_j.end(), vvms_index[v]) != Bubbles_p_i_j.end()) {
                            // Add random perturbation for Bubble case with a random seed
                            vvm::Init::RandomPerturbation(*local_vvms[local_v], local_vvms[local_v]->step + seed, -0.001, 0.001, 1.);
                        } else {
                            // Add random perturbation for Not Bubble case with the same random seed
                            vvm::Init::RandomPerturbation(*local_vvms[local_v], local_vvms[local_v]->step);
                        }
                    }
                    vvm::AddForcing(*local_vvms[local_v]);
                #endif

                vvm::BoundaryProcess2D_all(*local_vvms[local_v]);

                vvm::PoissonSolver::pubarTop_pt(*local_vvms[local_v]);
                vvm::PoissonSolver::cal_w(*local_vvms[local_v], p, i, j);
                vvm::PoissonSolver::cal_u(*local_vvms[local_v]);

                vvm::Iteration::updateMean(*local_vvms[local_v]);

                #if defined(DIFFUSION_VVM)
                    if (time_vvm == 0) std::cout << "Constant Diffusion" << std::endl;
                    vvm::NumericalProcess::DiffusionAll(*local_vvms[local_v]);
                #else
                    vvm::Turbulence::RKM_RKH(*local_vvms[local_v]);
                #endif

                vvm::NumericalProcess::Nudge_theta(*local_vvms[local_v]);
                if (local_vvms[local_v]->CASE != 2) vvm::NumericalProcess::Nudge_zeta(*local_vvms[local_v]);
                if (vvm_moisture_nudge_time != 0 && local_vvms[local_v]->CASE == 1) vvm::NumericalProcess::Nudge_qv(*local_vvms[local_v]);

                #if defined(WATER)
                    #if defined(KESSLER_MICROPHY)
                        vvm::MicroPhysics::autoconversion(*local_vvms[local_v]);
                        vvm::MicroPhysics::accretion(*local_vvms[local_v]);
                        vvm::MicroPhysics::evaporation(*local_vvms[local_v]);
                        vvm::MicroPhysics::condensation(*local_vvms[local_v]);
                        vvm::MicroPhysics::NegativeValueProcess(local_vvms[local_v]->qvp, local_vvms[local_v]->nx, local_vvms[local_v]->nz);
                        vvm::MicroPhysics::NegativeValueProcess(local_vvms[local_v]->qcp, local_vvms[local_v]->nx, local_vvms[local_v]->nz);
                        vvm::MicroPhysics::NegativeValueProcess(local_vvms[local_v]->qrp, local_vvms[local_v]->nx, local_vvms[local_v]->nz);
                    #endif

                    #if defined(P3_MICROPHY)

                    for (int i_vvm = 1; i_vvm < vvm_nx-1; i_vvm++) {
                        double *qcp1d = local_vvms[local_v]->qcp[i_vvm];
                        double *ncp1d = local_vvms[local_v]->ncp[i_vvm];
                        double *qrp1d = local_vvms[local_v]->qrp[i_vvm];
                        double *nrp1d = local_vvms[local_v]->nrp[i_vvm];
                        double *th1d = local_vvms[local_v]->th[i_vvm];
                        double *thp1d = local_vvms[local_v]->thp[i_vvm];
                        double *qv1d = local_vvms[local_v]->qv[i_vvm];
                        double *qvp1d = local_vvms[local_v]->qvp[i_vvm];
                        double *qitotp1d = local_vvms[local_v]->qitotp[i_vvm];
                        double *qirimp1d = local_vvms[local_v]->qirimp[i_vvm];
                        double *qiliqp1d = local_vvms[local_v]->qiliqp[i_vvm];
                        double *nip1d = local_vvms[local_v]->nip[i_vvm];
                        double *birimp1d = local_vvms[local_v]->birimp[i_vvm];
                        double *zi_all1d = local_vvms[local_v]->zi_all[i_vvm];
                        double *ssat_all1d = local_vvms[local_v]->ssat_all[i_vvm];
                        double *w_all1d = local_vvms[local_v]->w_all[i_vvm];
                        double *pb_all1d = local_vvms[local_v]->pb_all[i_vvm];
                        double *dz_all1d = local_vvms[local_v]->dz_all[i_vvm];
                        double *precip_liq1d = &local_vvms[local_v]->precip_liq[i_vvm];
                        double *precip_sol1d = &local_vvms[local_v]->precip_sol[i_vvm];
                        double *diag_ze1d = local_vvms[local_v]->diag_ze[i_vvm];
                        double *diag_effc1d = local_vvms[local_v]->diag_effc[i_vvm];
                        double *diag_effi1d = local_vvms[local_v]->diag_effi[i_vvm];
                        double *diag_vmi1d = local_vvms[local_v]->diag_vmi[i_vvm];
                        double *diag_di1d = local_vvms[local_v]->diag_di[i_vvm];
                        double *diag_rhoi1d = local_vvms[local_v]->diag_rhoi[i_vvm];
                        double *diag_2d1d = local_vvms[local_v]->diag_2d[i_vvm];
                        double *diag_3d1d = local_vvms[local_v]->diag_3d[i_vvm][0];
                        double *cldfrac1d = local_vvms[local_v]->cldfrac[i_vvm];

                        __microphy_p3_MOD_p3_main(
                            qcp1d, ncp1d, qrp1d, nrp1d, 
                            th1d, thp1d, qv1d, qvp1d, &local_vvms[local_v]->dt,
                            qitotp1d, qirimp1d, qiliqp1d, nip1d,
                            birimp1d, zi_all1d, ssat_all1d, w_all1d, pb_all1d,
                            dz_all1d, &local_vvms[local_v]->step, precip_liq1d, precip_sol1d, &one, &one, &one, &vvm_nz, 
                            &vvm::P3::nCat, diag_ze1d, diag_effc1d, diag_effi1d,
                            diag_vmi1d, diag_di1d, diag_rhoi1d, 
                            &vvm::P3::n_diag_2d, diag_2d1d, &vvm::P3::n_diag_3d, diag_3d1d,
                            &vvm::P3::log_predictNc, vvm::P3::model_name, &vvm::P3::clbfact_dep, 
                            &vvm::P3::clbfact_sub, &vvm::P3::debug_on, &vvm::P3::scpf_on, 
                            &vvm::P3::scpf_pfrac, &vvm::P3::scpf_resfact, cldfrac1d, 
                            &vvm::P3::trplMomI, &vvm::P3::liqfrac, 
                            nullptr, nullptr, nullptr, nullptr, nullptr,
                            nullptr, nullptr, nullptr, nullptr, nullptr,
                            nullptr, nullptr, nullptr, nullptr, nullptr, strlen(vvm::P3::model_name)
                        );
                        local_vvms[local_v]->precip[i_vvm] = local_vvms[local_v]->precip_sol[i_vvm] + local_vvms[local_v]->precip_liq[i_vvm];
                    }
                    #endif
                #endif

                vvm::BoundaryProcess2D_all(*local_vvms[local_v]);

                #if defined(TIMEFILTER) && !defined(AB2)
                    vvm::NumericalProcess::timeFilterAll(*local_vvms[local_v]);
                #endif

                // Apply large-scale forcing from master
                double total_heating = received_data[local_v].value / exchange_coeff;  // Use hp from scattered data
                for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
                    double heating = total_heating * heating_weight[k_vvm];
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        local_vvms[local_v]->thp[i_vvm][k_vvm] += coupling_vvm_param * heating * local_vvms[local_v]->dt;
                    }
                }

                // VVM next step
                vvm::Iteration::nextTimeStep(*local_vvms[local_v]);
                local_vvms[local_v]->step++;
            }
            time_vvm = local_vvms[0]->step * local_vvms[0]->dt;
        }

        // Gather th_mean from workers after VVM steps
        data_send* send_data = new data_send[local_size];
        for (int v = start; v <= end; v++) {
            int local_v = v - start;
            int p = vvms_index[v].p;
            int i = vvms_index[v].i;
            int j = vvms_index[v].j;

            double th_mean = 0.;
            for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
                for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                    th_mean += local_vvms[local_v]->th[i_vvm][k_vvm];
                }
            }
            th_mean /= ((vvm_nx-2) * k_couple);
            send_data[local_v] = {p, i, j, th_mean};
        }

        data_send* gathered_data = (rank == MASTER_RANK) ? new data_send[total_vvm_size] : nullptr;
        MPI_Gatherv(send_data, local_size, data_send_type,
                    gathered_data, gather_counts, displs_gather, data_send_type,
                    MASTER_RANK, MPI_COMM_WORLD);

        // Update th_mean_all on master with latest VVM state
        if (rank == MASTER_RANK) {
            for (int i = 0; i < total_vvm_size; i++) {
                int p = gathered_data[i].p;
                int x = gathered_data[i].i;
                int y = gathered_data[i].j;
                th_mean_all[p][x][y] = gathered_data[i].value;
            }
        }

        // Update CSSWM with new th_mean (master only)
        if (rank == MASTER_RANK) {
            for (int i = 0; i < total_vvm_size; i++) {
                int p = vvms_index[i].p;
                int x = vvms_index[i].i;
                int y = vvms_index[i].j;
                model_csswm->csswm[p].h[x][y] = th_mean_all[p][x][y] * exchange_coeff;
            }
        }

        next_coupling_time += Couple_time;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Cleanup
    for (int v = 0; v < local_size; v++) {
        delete local_vvms[v];
        delete config_vvms[v];
    }
    delete[] local_vvms;
    delete[] config_vvms;
    delete[] gather_counts;
    delete[] displs_gather;
    delete[] scatter_counts;
    delete[] displs_scatter;
    delete[] scattered_data;
    delete[] received_data;
    if (rank == MASTER_RANK) {
        delete[] th_mean_all;
        delete[] Q_all;
        delete[] q_all;
    }
    MPI_Type_free(&data_send_type);
    MPI_Finalize();
    return 0;
}

void output_forcing(std::string dir, int n, double q[6][NX][NY]) {
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

    NcVar q_all = dataFile.addVar("forcing", ncDouble, xyDim);

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
