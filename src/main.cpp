#include "allocate_csswm_vvms.hpp"
#include "reading_config.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif
#if defined(NCOUTPUT)
    #include <netcdf>
#endif
#include <cstdio>

using namespace netCDF;

void output_forcing(std::string dir, int n, double q[6][NX][NY]);

CSSWM model_csswm;

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

int main(int argc, char **argv) {
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

    // Read configuration file
    std::map<std::string, std::string> config = read_config("../config.txt");
    std::cout << "Configuration key-value pairs:" << std::endl;
    for (const auto& pair : config) {
        std::cout << pair.first << " = " << pair.second << std::endl;
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

    printf("vvm_moisutre_nudge_time: %f\n", vvm_moisture_nudge_time);

    model_csswm.output_path = path + "csswm/";
    model_csswm.gravity = csswm_gravity;
    model_csswm.dt = csswm_dt;
    model_csswm.timeend = csswm_timeend;
    model_csswm.outputstep = csswm_outputstep;
    model_csswm.diffusion_kx = csswm_diffusion_kx;
    model_csswm.diffusion_ky = csswm_diffusion_ky;
    model_csswm.diffusion_ts = csswm_diffusion_ts;
    model_csswm.addforcingtime = csswm_addforcingtime;
    model_csswm.csswm_h_nudge_time = csswm_h_nudge_time;
    
    printf("output path: %s\n", path.c_str());
    printf("seed: %d\n", seed);
    printf("Coupling time: %f\n", Couple_time);
    printf("Gravity wave speed: %f\n", model_csswm.gravity);
    // Display the parsed tuples
    for (const vvm_index& tuple : Bubbles_p_i_j) {
        printf("Bubble: p=%d, i=%d, j=%d\n", tuple.p, tuple.i, tuple.j);
    }
    for (const vvm_index& tuple : NotBubbles_p_i_j) {
        printf("Not Bubble: p=%d, i=%d, j=%d\n", tuple.p, tuple.i, tuple.j);
    }
    printf("vvm_dt: %f s\n", vvm_dt);
    printf("vvm_dx: %f m\n", vvm_dx);
    printf("vvm_dz: %f m\n", vvm_dz);

    CSSWM::Init::Init2d(model_csswm);

    Config_VVM**** config_vvms = allocate_and_initialize_config(6, NX, NY);
    std::string path_vvm;

    // Make configurations for all VVMs
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                path_vvm = path + "vvm/" + std::to_string(p) + "_" + std::to_string(i) + "_" + std::to_string(j) + "/";
                config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 0, vvm_xrange, vvm_zrange, vvm_dx, vvm_dz, vvm_dt, vvm_timeend, vvm_outputstep, vvm_moisture_nudge_time));
            }
        }
    }
    // Change some configurations for Bubbles
    for (auto Bubble : Bubbles_p_i_j) {
        path_vvm = path + "vvm/" + std::to_string(Bubble.p) + "_" + std::to_string(Bubble.i) + "_" + std::to_string(Bubble.j) + "/";
        config_vvms[Bubble.p][Bubble.i][Bubble.j] =  new Config_VVM(createConfig(path_vvm, -1, Bubble_case, vvm_xrange, vvm_zrange, vvm_dx, vvm_dz, vvm_dt, vvm_timeend, vvm_outputstep, vvm_moisture_nudge_time));
    }
    printf("Configurations are set.\n");

    int total_size = Bubbles_p_i_j.size() + NotBubbles_p_i_j.size();
    vvm_index vvms_index[total_size];
    int count = 0;
    for (int size = 0; size < (int) Bubbles_p_i_j.size(); size++) {
        vvms_index[count] = Bubbles_p_i_j[size];
        count++;
    }
    for (int size = 0; size < (int) NotBubbles_p_i_j.size(); size++) {
        vvms_index[count] = NotBubbles_p_i_j[size];
        count++;
    }
    printf("count: %d\n", count);
    if (count != total_size) {
        printf("Error: count != total_size\n");
        return 1;
    }

    for (int size = 0; size < total_size; size++) {
        printf("p: %d, i: %d, j: %d\n", vvms_index[size].p, vvms_index[size].i, vvms_index[size].j);
    }
    
    vvm**** vvms = allocate_and_initialize(6, NX, NY);
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

    double temp_csswm = csswm_timeend / csswm_dt, temp_vvm = vvm_timeend / vvm_dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    CSSWM::Outputs::create_all_directory(model_csswm);
    // create Q_all directory
    CSSWM::Outputs::create_directory(model_csswm.output_path + (std::string) "Q_all/");
    CSSWM::Outputs::create_directory(model_csswm.output_path + (std::string) "q_all/");

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

    // Copy grads ctl file to the output directory
    std::string src1 = "../CSSWM/scripts/csswm.ctl";
    std::string des1 = model_csswm.output_path + "nc/csswm.ctl";

    // Construct the command
    std::string cmd1 = "cp " + src1 + " " + des1;

    // Execute the command
    system(cmd1.c_str());

    std::string src2 = "../2DVVM/scripts/vvm.ctl";
    std::string des2;
    std::string cmd2;

    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;
        des2 = vvms[p][i][j]->outputpath + "nc/vvm.ctl";
        cmd2 = "cp " + src2 + " " + des2;
        system(cmd2.c_str());
    }


    int vvm_nx = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->nx;
    int vvm_nz = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->nz;

    double exchange_coeff = 0.;
    double coupling_vvm_param = 1.;

    double th_mean = 0.;
    double th_mean_all[6][NX][NY];
    double Q_all[6][NX][NY];
    double q_all[6][NX][NY];

    int k_couple = vvm_nx - 2;
    k_couple = 12000. / vvms[NotBubbles_p_i_j[0].p][NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j]->dz;

    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;

        printf("p: %d, i: %d, j: %d\n", p, i, j);
        
        th_mean_all[p][i][j] = 0.;
        for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                th_mean_all[p][i][j] += vvms[p][i][j]->th[i_vvm][k_vvm];
            }
        }
        th_mean_all[p][i][j] /= ((vvm_nx-2) * k_couple);
    }
    exchange_coeff = model_csswm.csswm[NotBubbles_p_i_j[0].p].h[NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j] / th_mean_all[NotBubbles_p_i_j[0].p][NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j];

    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;

        double total_heating = (model_csswm.csswm[p].h[i][j] / exchange_coeff - th_mean_all[p][i][j]) * (vvm_nz-2);

        for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
            double heating = total_heating * heating_weight[k_vvm];
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                vvms[p][i][j]->th[i_vvm][k_vvm] += heating;
                vvms[p][i][j]->thm[i_vvm][k_vvm] = vvms[p][i][j]->th[i_vvm][k_vvm];
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

    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;
        for (int k_vvm = 0; k_vvm <= vvm_nz-1; k_vvm++) {
            for (int i_vvm = 0; i_vvm <= vvm_nx-1; i_vvm++) {
                vvms[p][i][j]->dz_all[i_vvm][k_vvm] = vvms[p][i][j]->dz;
                vvms[p][i][j]->w_all[i_vvm][k_vvm] = 0.;
                vvms[p][i][j]->pb_all[i_vvm][k_vvm] = vvms[p][i][j]->pb[k_vvm];
                vvms[p][i][j]->zi_all[i_vvm][k_vvm] = 0.;
                vvms[p][i][j]->ssat_all[i_vvm][k_vvm] = 0.;
            }
        }
        for (int k_vvm = 0; k_vvm < vvm_nz; k_vvm++) {
            for (int i_vvm = 0; i_vvm < vvm_nx; i_vvm++) {
                vvms[p][i][j]->qiliqp[i_vvm][k_vvm] = 0.;
            }
        }
    }
    int one = 1;
    #endif
    // initialize Q_all, q_all
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                th_mean_all[p][i][j] = 0.;
                Q_all[p][i][j] = 0.;
                q_all[p][i][j] = 0.;
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif

    double next_coupling_time = Couple_time;
    while (vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step < nmax_vvm && model_csswm.step < nmax_csswm) {

        double time_vvm = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step * vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->dt;
        double time_csswm = model_csswm.step * DT;

        if (model_csswm.csswm[NotBubbles_p_i_j[0].p].h[NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j] != model_csswm.csswm[NotBubbles_p_i_j[0].p].h[NotBubbles_p_i_j[0].i][NotBubbles_p_i_j[0].j]) {
            printf("Nan\n");
            return 1;
        }

        while (next_coupling_time != time_csswm) {
            printf("csswm_step: %d, csswm_time: %f\n", model_csswm.step, time_csswm);

            // Output for CSSWM
            if (model_csswm.step % model_csswm.outputstep == 0 || model_csswm.step == model_csswm.timeend-1 || model_csswm.step == model_csswm.timeend-2) {
                #ifdef NCOUTPUT
                    CSSWM::Outputs::huv_nc(model_csswm.step, model_csswm);
                #endif
            }

            #if defined(EquatorialWave)
                if (model_csswm.step * model_csswm.dt >= model_csswm.addforcingtime) model_csswm.status_add_forcing = false;
                else model_csswm.status_add_forcing = true;
            #endif

            // Prediction for CSSWM
            CSSWM::Iteration::ph_pt_4(model_csswm);
            CSSWM::Iteration::pu_pt_4(model_csswm);
            CSSWM::Iteration::pv_pt_4(model_csswm);

            CSSWM::BP_h(model_csswm);
            model_csswm.BP_wind_interpolation2(model_csswm);

            #if defined(DIFFUSION)
                CSSWM::NumericalProcess::DiffusionAll(model_csswm);
            #endif

            model_csswm.BP_h(model_csswm);
            model_csswm.BP_wind_interpolation2(model_csswm);
            
            #if defined(TIMEFILTER) && !defined(AB2Time)
                CSSWM::NumericalProcess::timeFilterAll(model_csswm);
            #endif

            if (model_csswm.csswm_h_nudge_time != 0) {
                CSSWM::NumericalProcess::NudgeH(model_csswm);
            }

            CSSWM::Iteration::nextTimeStep(model_csswm);
            model_csswm.step++;
            time_csswm = model_csswm.step * DT;
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        // Get th_mean at time step n, which is before the iteration for CRM
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;

            th_mean = 0.;
            for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
                for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                    th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                }
            }
            th_mean /= ((vvm_nx-2) * k_couple);
            th_mean_all[p][i][j] = th_mean;

            q_all[p][i][j] = (model_csswm.csswm[p].hp[i][j] / exchange_coeff - th_mean_all[p][i][j]) / Couple_time;
            // Q_all[p][i][j] = (exchange_coeff * th_mean_all[p][i][j] - model_csswm.csswm[p].h[i][j]) / Couple_time;
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        std::cout << th_mean_all[1][46][47] << ", " << th_mean_all[1][47][47] << std::endl;

        output_forcing(model_csswm.output_path + (std::string) "q_all/", (int) next_coupling_time / Couple_time, q_all);

        while (time_vvm < next_coupling_time) {
            printf("VVM step: %d, time: %f\n", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, time_vvm);
            
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int size = 0; size < total_size; size++) {
                int p = vvms_index[size].p;
                int i = vvms_index[size].i;
                int j = vvms_index[size].j;

                if (vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step % vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->OUTPUTSTEP == 0) {
                    #if defined(OUTPUTTXT)
                        vvm::Output::outputalltxt(vvms[p][i][j]->step, *vvms[p][i][j]);
                    #endif

                    #pragma omp critical
                    {
                        #if defined(OUTPUTNC)
                            vvm::Output::output_nc(vvms[p][i][j]->step, *vvms[p][i][j]);
                        #endif
                    }
                }

                vvm::Iteration::pzeta_pt(*vvms[p][i][j]);
                vvm::Iteration::pth_pt(*vvms[p][i][j]);
                #if defined(WATER)
                    #if defined(KESSLER_MICROPHY)
                        vvm::Iteration::pqv_pt(*vvms[p][i][j]);
                        vvm::Iteration::pqc_pt(*vvms[p][i][j]);
                        vvm::Iteration::pqr_pt(*vvms[p][i][j]);
                    #endif

                    #if defined(P3_MICROPHY)
                        vvm::Iteration::pqmicrophy_pt(*vvms[p][i][j]);
                    #endif

                    if (vvms[p][i][j]->step * vvms[p][i][j]->dt <= vvms[p][i][j]->addforcingtime) vvms[p][i][j]->status_for_adding_forcing = true;
                    else vvms[p][i][j]->status_for_adding_forcing = false;

                    // Generate new random th perturbation for tropical forcing case
                    if (vvms[p][i][j]->status_for_adding_forcing == true) {
                        if (is_value_in_vvm_index(Bubbles_p_i_j, p, i, j)) {
                            // Add random perturbation for Bubble case with a random seed
                            vvm::Init::RandomPerturbation(*vvms[p][i][j], vvms[p][i][j]->step+seed, -0.001, 0.001, 1.);
                        }
                        else {
                            // Add random perturbation for Not Bubble case with the same random seed
                            vvm::Init::RandomPerturbation(*vvms[p][i][j], vvms[p][i][j]->step);
                        }
                    }
                    vvm::AddForcing(*vvms[p][i][j]);
                #endif
                vvm::BoundaryProcess2D_all(*vvms[p][i][j]);

                vvm::PoissonSolver::pubarTop_pt(*vvms[p][i][j]);
                vvm::PoissonSolver::cal_w(*vvms[p][i][j], p, i, j);
                vvm::PoissonSolver::cal_u(*vvms[p][i][j]);
                
                vvm::Iteration::updateMean(*vvms[p][i][j]);
                #if defined(DIFFUSION_VVM)
                    if (time_vvm == 0) std::cout << "Constant Diffusion" << std::endl;
                    vvm::NumericalProcess::DiffusionAll(*vvms[p][i][j]);
                #else
                    vvm::Turbulence::RKM_RKH(*vvms[p][i][j]);
                #endif
                vvm::NumericalProcess::Nudge_theta(*vvms[p][i][j]);
                if (vvms[p][i][j]->CASE != 2) vvm::NumericalProcess::Nudge_zeta(*vvms[p][i][j]);
                if (vvm_moisture_nudge_time != 0 && vvms[p][i][j]->CASE == 1) vvm::NumericalProcess::Nudge_qv(*vvms[p][i][j]);


                #if defined(WATER)
                    #if defined(KESSLER_MICROPHY)
                        vvm::MicroPhysics::autoconversion(*vvms[p][i][j]);
                        vvm::MicroPhysics::accretion(*vvms[p][i][j]);
                        vvm::MicroPhysics::evaporation(*vvms[p][i][j]);
                        vvm::MicroPhysics::condensation(*vvms[p][i][j]); // saturation adjustment

                        // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
                        vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qvp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
                        vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qcp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
                        vvm::MicroPhysics::NegativeValueProcess(vvms[p][i][j]->qrp, vvms[p][i][j]->nx, vvms[p][i][j]->nz);
                    #endif

                    #if defined(P3_MICROPHY)

                    for (int i_vvm = 1; i_vvm < vvm_nx-1; i_vvm++) {
                        double *qcp1d = vvms[p][i][j]->qcp[i_vvm];
                        double *ncp1d = vvms[p][i][j]->ncp[i_vvm];
                        double *qrp1d = vvms[p][i][j]->qrp[i_vvm];
                        double *nrp1d = vvms[p][i][j]->nrp[i_vvm];
                        double *th1d = vvms[p][i][j]->th[i_vvm];
                        double *thp1d = vvms[p][i][j]->thp[i_vvm];
                        double *qv1d = vvms[p][i][j]->qv[i_vvm];
                        double *qvp1d = vvms[p][i][j]->qvp[i_vvm];
                        double *qitotp1d = vvms[p][i][j]->qitotp[i_vvm];
                        double *qirimp1d = vvms[p][i][j]->qirimp[i_vvm];
                        double *qiliqp1d = vvms[p][i][j]->qiliqp[i_vvm];
                        double *nip1d = vvms[p][i][j]->nip[i_vvm];
                        double *birimp1d = vvms[p][i][j]->birimp[i_vvm];
                        double *zi_all1d = vvms[p][i][j]->zi_all[i_vvm];
                        double *ssat_all1d = vvms[p][i][j]->ssat_all[i_vvm];
                        double *w_all1d = vvms[p][i][j]->w_all[i_vvm];
                        double *pb_all1d = vvms[p][i][j]->pb_all[i_vvm];
                        double *dz_all1d = vvms[p][i][j]->dz_all[i_vvm];
                        double *precip_liq1d = &vvms[p][i][j]->precip_liq[i_vvm];
                        double *precip_sol1d = &vvms[p][i][j]->precip_sol[i_vvm];
                        double *diag_ze1d = vvms[p][i][j]->diag_ze[i_vvm];
                        double *diag_effc1d = vvms[p][i][j]->diag_effc[i_vvm];
                        double *diag_effi1d = vvms[p][i][j]->diag_effi[i_vvm];
                        double *diag_vmi1d = vvms[p][i][j]->diag_vmi[i_vvm];
                        double *diag_di1d = vvms[p][i][j]->diag_di[i_vvm];
                        double *diag_rhoi1d = vvms[p][i][j]->diag_rhoi[i_vvm];
                        double *diag_2d1d = vvms[p][i][j]->diag_2d[i_vvm];
                        double *diag_3d1d = vvms[p][i][j]->diag_3d[i_vvm][0];
                        double *cldfrac1d = vvms[p][i][j]->cldfrac[i_vvm];

                        __microphy_p3_MOD_p3_main(
                            qcp1d, ncp1d, qrp1d, nrp1d, 
                            th1d, thp1d, qv1d, qvp1d, &vvms[p][i][j]->dt,
                            qitotp1d, qirimp1d, qiliqp1d, nip1d,
                            birimp1d, zi_all1d, ssat_all1d, w_all1d, pb_all1d,
                            dz_all1d, &vvms[p][i][j]->step, precip_liq1d, precip_sol1d, &one, &one, &one, &vvm_nz, 
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
                        vvms[p][i][j]->precip[i_vvm] = vvms[p][i][j]->precip_sol[i_vvm] + vvms[p][i][j]->precip_liq[i_vvm];
                    }
                    #endif
                #endif

                vvm::BoundaryProcess2D_all(*vvms[p][i][j]);

                #if defined(TIMEFILTER) && !defined(AB2)
                    vvm::NumericalProcess::timeFilterAll(*vvms[p][i][j]);
                #endif

                // Large Scale Forcing
                double total_heating = q_all[p][i][j] * k_couple;
                double heating = 0.;

                for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
                    heating = total_heating * heating_weight[k_vvm];
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * heating * vvms[p][i][j]->dt;
                    }
                }

                // VVM next step
                vvm::Iteration::nextTimeStep(*vvms[p][i][j]);
                vvms[p][i][j]->step++;
            }
            #ifdef _OPENMP
            #pragma omp barrier
            #endif
            time_vvm = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step * vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->dt;
        }

        // Next time step for CSSWM (the iteration is done after adding the small scale forcing)
        #ifdef _OPENMP
        #pragma omp parallel for 
        #endif
        for (int size = 0; size < total_size; size++) {
            int p = vvms_index[size].p;
            int i = vvms_index[size].i;
            int j = vvms_index[size].j;

            double th_mean = 0.;
            for (int k_vvm = 1; k_vvm <= k_couple; k_vvm++) {
                for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                    th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                }
            }
            th_mean /= ((vvm_nx-2) * k_couple);

            Q_all[p][i][j] = (th_mean * exchange_coeff - model_csswm.csswm[p].h[i][j]);
            model_csswm.csswm[p].h[i][j] = th_mean * exchange_coeff;
            
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        output_forcing(model_csswm.output_path + (std::string) "Q_all/", (int) next_coupling_time / Couple_time, Q_all);

        next_coupling_time += Couple_time;
    }

    deallocate_config(config_vvms, 6, NX, NY);
    deallocate(vvms, 6, NX, NY);
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
