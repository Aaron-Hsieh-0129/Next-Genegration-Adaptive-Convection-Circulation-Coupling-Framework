#include "../2DVVM/src/Declare.hpp"
#include "../CSSWM/src/construction.hpp"
#include "../CSSWM/src/define.hpp"
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <netcdf>

using namespace netCDF;

#define PROFILE
// #define AB2_Couple
#define Couple_10km


// CASE0: Nothing, CASE1:Bubble
Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double Kx, double Kz) {
    return Config_VVM(3.0, 200.0, 200.0, 100000, 20000, 1500000.0, 
                      10000, path, 50, 
                      Kx, Kz, 0.01, 1E-22, 9.80665, 1003.5, 716.5, 287.0, 
                      2.5E6, 1E5, 96500.0, addforcingtime, CASE);
}



// void output_qall(std::string dir,int n, double q[6][NX][NY]);
// void output_Qall(std::string dir,int n, double Q[6][NX][NY]);

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


CSSWM model_csswm;

int main(int argc, char **argv) {
    #if defined(PROFILE)
        // This heating weight follows the Q1 heating profile for the data in 2DVVM/input/init.txt
        #if defined(Couple_10km)
        // double heating_weight[52] = {
        //     0.,
        //     0.008147876891322802, 0.018476288762754464, 0.0266324883066731, 0.032458345123757844, 0.035953859214008684, 0.03878356109659271, 0.04061454466767648, 0.0423623017128019, 0.0436106996021772, 0.04485909749155251, 0.04544168317326098, 
        //     0.0459410423290111, 0.045774589277094394, 0.045275230121344275, 0.04477587096559415, 0.0436106996021772, 0.04261198129067696, 0.04136358340130166, 0.03969905288213459, 0.03811774888892588, 0.03578740616209198, 
        //     0.03287447775354961, 0.02996154934500724, 0.026466035254756395, 0.022637615060672137, 0.01889242139254623, 0.015646586880170447, 0.012733658471628075, 0.009903956589044059, 0.008089618323151955, 0.006516636982539074, 
        //     0.0049686235997137, 0.0039282920252342824, 0.002904605755946535, 0.001989113970404647, 0.0013898829835045025, 0.000798974649200193, 0.0, 0.0, 0.0, 0.0, 
        //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        //     0.
        // };

        double heating_weight[52] = {
            0.0, 
            0.0013851364579662475, 0.004791691139042248, 0.008173773444535806, 0.010865734758604487, 0.01350875132150828, 0.015662320372763223, 0.017375386663534203, 0.01908845295430518, 0.020214182231097534, 0.02114413250323035, 0.022025138024198282, 
            0.022808254042836443, 0.02334664630565018, 0.023885038568463916, 0.02442343083127765, 0.0249128783429265, 0.025255491601080697, 0.025647049610399777, 0.026038607619718858, 0.026381220877873056, 0.026576999882532593, 
            0.026723834136027247, 0.026870668389521905, 0.02701750264301656, 0.027115392145346327, 0.02691961314068679, 0.026772778887192136, 0.02662594463369748, 0.026479110380202824, 0.02633227612670817, 0.025989662868553975, 
            0.025647049610399777, 0.025353381103410465, 0.025059712596421157, 0.024717099338266962, 0.024325541328947882, 0.02383609381729903, 0.02334664630565018, 0.022857198794001325, 0.022416696033517362, 0.020378814575924876, 
            0.01834093311833239, 0.016303051660739903, 0.014265170203147413, 0.012227288745554924, 0.01018940728796244, 0.008151525830369951, 0.006113644372777461, 0.004075762915184976, 0.002037881457592486,
            0.0
        };
        #else
        double heating_weight[102] = {
            0.,
            0.001177598498000163, 0.004073741800502331, 0.00694907947583135, 0.009237698464877604, 0.011484706199577562, 0.013315601390814563, 0.014771995292934909, 0.01622838919505525, 0.017185448045020046, 0.01797606187759966, 0.018725064455832982, 
            0.01939084452537371, 0.01984856832318296, 0.02030629212099221, 0.020764015918801462, 0.021180128462264414, 0.021471407242688485, 0.021804297277458848, 0.02213718731222921, 0.022428466092653282, 0.022594911110038463, 
            0.02271974487307735, 0.022844578636116237, 0.022969412399155124, 0.023052634907847713, 0.02288618989046253, 0.022761356127423645, 0.022636522364384758, 0.02251168860134587, 0.022386854838306984, 0.02209557605788292, 
            0.021804297277458848, 0.021554629751381074, 0.0213049622253033, 0.021013683444879232, 0.02068079341010887, 0.020264680866645915, 0.01984856832318296, 0.01943245577972, 0.019057954490603345, 0.01864184194714039, 
            0.017892839368907072, 0.01714383679067375, 0.016436445466786725, 0.015687442888553407, 0.014980051564666384, 0.01418943773208677, 0.013232378882121974, 0.012275320032157176, 0.011318261182192379, 0.010361202332227582, 
            0.009445754736609082, 0.008571918395336876, 0.007822915817103556, 0.0071155244932165325, 0.006366521914983213, 0.005659130591096189, 0.0049517392672091655, 0.004452404215053619, 0.004044613922459924, 0.0036493070061701166, 
            0.0032581612153149385, 0.002867015424459761, 0.002484191884473842, 0.0022220409820921804, 0.001964051205145148, 0.0017060614281981159, 0.0014522327766857135, 0.0011984041251733107, 0.0009945089788764628, 0.0008405473377951694, 
            0.0006949079475831351, 0.0005451074319364712, 0.0003994680417244369, 0.00025632532677318036, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.
        };
        #endif
    #endif

    omp_set_num_threads(128);
    Eigen::setNbThreads(1);

    std::string path = "/data/Aaron/TMIF/0901_stable_version/dt600_1_csswm_1_vvm_2E5diff_7vvm_3B_4non_10kmcouple/";

    // Config_CSSWM config_csswm(1800., 1., 1., 0.1, 86400 * 3 * 24., path + "csswm/", 
    //                     1, 2E5, 2E5, 0.06, 1200. * 60.);
    // CSSWM model_csswm(config_csswm);


    
    CSSWM::Init::Init2d(model_csswm);
    

    Config_VVM**** config_vvms = allocate_and_initialize_config(6, NX, NY);
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
                // if (p == 1 && (NY/2-10 <= j && j <= NY/2+10) && (i == NX/2)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));
                // if (p == 1 && (j == NY/2) && (38 <= i && i <= 56)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));
                // if (p == 1 && (j == NY/2) && (34 <= i && i <= 61)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));
                // else config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 0, 70, 70));

                // if (p == 1 && i == 46 && (38 <= j && j <= 56 && j != 47)) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));
                // else config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 0, 70, 70));

                // if (p == 1 && i == NX/2 && j == NY/2) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));

                if (p == 1 && (46 <= i && i <= 48) && j == 47) config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, -1, 1, 200, 200));
                else config_vvms[p][i][j] = new Config_VVM(createConfig(path_vvm, 10, 0, 70, 70));
            }
        }
    }
    printf("Configurations are set.\n");

    int total_size = 7;
    vvm_index vvms_index[total_size];
    int count = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i <= NX-2; i++) {
            for (int j = 2; j <= NY-2; j++) {
                if (p == 1 && (44 <= i && i <= 50) && j == 47) {
                    vvms_index[count] = {p, i, j};
                    count++;
                }
            }
        }
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

    double temp_csswm = TIMEEND / DT, temp_vvm = TIMEEND / config_vvms[1][NX/2][NY/2]->dt;
    int nmax_csswm = (int) temp_csswm, nmax_vvm = (int) temp_vvm;

    CSSWM::Outputs::create_all_directory();
    // create Q_all directory
    CSSWM::Outputs::create_directory(OUTPUTPATH + (std::string) "Q_all/");

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
    double th_mean_all[6][NX][NY];
    #if defined(AB2_Couple)
        double Q_all[2][6][NX][NY];
        double q_all[2][6][NX][NY];
    #else
        double Q_all[6][NX][NY];
        double q_all[6][NX][NY];
    #endif

    for (int size = 0; size < total_size; size++) {
        int p = vvms_index[size].p;
        int i = vvms_index[size].i;
        int j = vvms_index[size].j;

        printf("p: %d, i: %d, j: %d\n", p, i, j);
        
        th_mean_all[p][i][j] = 0.;
        for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                th_mean_all[p][i][j] += vvms[p][i][j]->th[i_vvm][k_vvm];
            }
        }
        th_mean_all[p][i][j] /= ((vvm_nx-2) * (vvm_nz-2));
        

        double total_heating = (model_csswm.csswm[p].h[i][j] / exchange_coeff - th_mean_all[p][i][j]) * (vvm_nz-2);

        for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
            double heating = total_heating * heating_weight[k_vvm];
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                vvms[p][i][j]->th[i_vvm][k_vvm] += heating;
                vvms[p][i][j]->thm[i_vvm][k_vvm] = vvms[p][i][j]->th[i_vvm][k_vvm];
            }
        }   
    }

    double tmp_th = 0.;
    for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
        for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
            tmp_th += vvms[1][47][47]->th[i_vvm][k_vvm];
        }
    }
    tmp_th /= ((vvm_nx-2) * (vvm_nz-2));

    printf("tmp_th: %f, psuedo_h: %f\n", tmp_th, tmp_th*exchange_coeff);
    printf("h: %f\n", model_csswm.csswm[1].h[47][47]);

    #if defined(Couple_10km)
        double tmp2_th = 0.;
        int k_10km = 10000. / vvms[1][47][47]->dz;
        printf("z: %f\n", vvms[1][47][47]->z[k_10km]);
        for (int k_vvm = 1; k_vvm <= k_10km; k_vvm++) {
            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                tmp2_th += vvms[1][47][47]->th[i_vvm][k_vvm];
            }
        }
        tmp2_th /= ((vvm_nx-2) * k_10km);
        double exchange_coeff2 = exchange_coeff * tmp_th / tmp2_th;
    #endif


    // initialize Q_all, q_all
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                th_mean_all[p][i][j] = 0.;
                #if defined(AB2_Couple)
                    Q_all[0][p][i][j] = Q_all[1][p][i][j] = 0.;
                    q_all[0][p][i][j] = q_all[1][p][i][j] = 0.;
                #else
                    Q_all[p][i][j] = 0.;
                    q_all[p][i][j] = 0.;
                #endif
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif

    while (vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step < nmax_vvm || model_csswm.step < nmax_csswm) {

        double time_vvm = vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step * vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->dt;
        double time_csswm = model_csswm.step * DT;
        printf("n_vvm: %d, time_vvm: %f, n_csswm: %d,  time_csswm: %f\n", vvms[vvms_index[0].p][vvms_index[0].i][vvms_index[0].j]->step, time_vvm, model_csswm.step, time_csswm);

        if (model_csswm.csswm[1].h[47][47] != model_csswm.csswm[1].h[47][47]) {
            printf("Nan\n");
            return 1;
        }
        if (time_vvm == time_csswm) {
            // Exchange information for small scale forcing
            if (time_csswm != 0) {
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for (int size = 0; size < total_size; size++) {
                    int p = vvms_index[size].p;
                    int i = vvms_index[size].i;
                    int j = vvms_index[size].j;

                    double th_mean = 0.;
                    #if defined(Couple_10km)
                        for (int k_vvm = 1; k_vvm <= k_10km; k_vvm++) {
                            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                                th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                            }
                        }
                        th_mean /= ((vvm_nx-2) * k_10km);
                        model_csswm.csswm[p].hp[i][j] = th_mean * exchange_coeff2;
                    #else
                        for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                            for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                                th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                            }
                        }
                        th_mean /= ((vvm_nx-2) * (vvm_nz-2));
                        // model_csswm.csswm[p].hp[i][j] += coupling_csswm_param * Q_all[p][i][j] * DT;
                        model_csswm.csswm[p].hp[i][j] = th_mean * exchange_coeff;
                    #endif
                    
                    #if defined(AB2_Couple)
                        model_csswm.csswm[p].hp[i][j] += coupling_csswm_param * (1.5*Q_all[(model_csswm.step+1)%2][p][i][j] - 0.5*Q_all[model_csswm.step%2][p][i][j]) * DT;
                    #endif
                }
                #ifdef _OPENMP
                #pragma omp barrier
                #endif

                // After adding the small scale forcing, we need to iterate the CSSWM model
                CSSWM::Iteration::nextTimeStep(model_csswm);
            }

            // Output for CSSWM
            if (model_csswm.step % OUTPUTINTERVAL == 0 || model_csswm.step == TIMEEND-1 || model_csswm.step == TIMEEND-2) {
                #ifdef NCOUTPUT
                    CSSWM::Outputs::huv_nc(model_csswm.step, model_csswm);
                #endif
            }

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
        }


        if (time_csswm == time_vvm) {
            // Get th_mean at time step n, which is before the iteration for CRM
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int size = 0; size < total_size; size++) {
                int p = vvms_index[size].p;
                int i = vvms_index[size].i;
                int j = vvms_index[size].j;
                #if defined(Couple_10km)
                    th_mean = 0.;
                    for (int k_vvm = 1; k_vvm <= k_10km; k_vvm++) {
                        for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                            th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                        }
                    }
                    th_mean /= ((vvm_nx-2) * k_10km);
                #else
                    th_mean = 0.;
                    for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                        for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                            th_mean += vvms[p][i][j]->th[i_vvm][k_vvm];
                        }
                    }
                    th_mean /= ((vvm_nx-2) * (vvm_nz-2));
                #endif
                
                th_mean_all[p][i][j] = th_mean;

                #if defined(AB2_Couple)
                    q_all[(model_csswm.step+1)%2][p][i][j] = (model_csswm.csswm[p].hp[i][j] / exchange_coeff - th_mean_all[p][i][j]) / DT;
                    Q_all[(model_csswm.step+1)%2][p][i][j] = (exchange_coeff * th_mean_all[p][i][j] - model_csswm.csswm[p].h[i][j]) / DT;
                #else
                    #if defined(Couple_10km)
                        q_all[p][i][j] = (model_csswm.csswm[p].hp[i][j] / exchange_coeff2 - th_mean_all[p][i][j]) / DT;
                        Q_all[p][i][j] = (exchange_coeff2 * th_mean_all[p][i][j] - model_csswm.csswm[p].h[i][j]) / DT;
                    #else
                        q_all[p][i][j] = (model_csswm.csswm[p].hp[i][j] / exchange_coeff - th_mean_all[p][i][j]) / DT;
                        Q_all[p][i][j] = (exchange_coeff * th_mean_all[p][i][j] - model_csswm.csswm[p].h[i][j]) / DT;
                        // Q_all[p][i][j] = (exchange_coeff * th_mean_all[p][i][j] - model_csswm.csswm[p].hp[i][j]) / DT;
                    #endif
                #endif
            }
            #ifdef _OPENMP
            #pragma omp barrier
            #endif
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
            vvm::NumericalProcess::Nudge_theta(*vvms[p][i][j]);
            vvm::NumericalProcess::Nudge_zeta(*vvms[p][i][j]);

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

            // Large Scale Forcing
            #if defined(PROFILE)
                #if defined(AB2_Couple)
                    double total_heating1 = q_all[model_csswm.step%2][p][i][j] * (vvm_nz-2);
                    double total_heating2 = q_all[(model_csswm.step+1)%2][p][i][j] * (vvm_nz-2);
                    double heating1 = 0.;
                    double heating2 = 0.;
                #else
                    #if defined(Couple_10km)
                        double total_heating = q_all[p][i][j] * k_10km;
                    #else
                        double total_heating = q_all[p][i][j] * (vvm_nz-2);
                    #endif
                    double heating = 0.;
                #endif
            #endif

            #if defined(Couple_10km)
                for (int k_vvm = 1; k_vvm <= k_10km; k_vvm++) {
                    #if defined(PROFILE)
                        #if defined(AB2_Couple)
                            heating1 = total_heating1 * heating_weight[k_vvm];
                            heating2 = total_heating2 * heating_weight[k_vvm];
                        #else
                            heating = total_heating * heating_weight[k_vvm];
                        #endif
                    #endif
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        #if defined(PROFILE)
                            #if defined(AB2_Couple)
                                vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * (1.5*heating2 - 0.5*heating1) * vvms[p][i][j]->dt;
                            #else
                                vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * heating * vvms[p][i][j]->dt;
                            #endif
                        #else
                            vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * vvms[p][i][j]->dt * q_all[p][i][j];
                        #endif
                    }
                }
            #else
                for (int k_vvm = 1; k_vvm <= vvm_nz-2; k_vvm++) {
                    #if defined(PROFILE)
                        #if defined(AB2_Couple)
                            heating1 = total_heating1 * heating_weight[k_vvm];
                            heating2 = total_heating2 * heating_weight[k_vvm];
                        #else
                            heating = total_heating * heating_weight[k_vvm];
                        #endif
                    #endif
                    for (int i_vvm = 1; i_vvm <= vvm_nx-2; i_vvm++) {
                        #if defined(PROFILE)
                            #if defined(AB2_Couple)
                                vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * (1.5*heating2 - 0.5*heating1) * vvms[p][i][j]->dt;
                            #else
                                vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * heating * vvms[p][i][j]->dt;
                            #endif
                        #else
                            vvms[p][i][j]->thp[i_vvm][k_vvm] += coupling_vvm_param * vvms[p][i][j]->dt * q_all[p][i][j];
                        #endif
                    }
                }
            #endif

            // VVM next step
            vvm::Iteration::nextTimeStep(*vvms[p][i][j]);
            vvms[p][i][j]->step++;
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        // Next time step for CSSWM (the iteration is done after adding the small scale forcing)
        if (time_csswm == time_vvm) {
            model_csswm.step++;
        }
    }

    deallocate_config(config_vvms, 6, NX, NY);
    deallocate(vvms, 6, NX, NY);
    return 0;
}

// void output_qall(std::string dir, int n, double q[6][NX][NY]) {
//     NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", NX);
//     NcDim yDim = dataFile.addDim("y", NY);
//     NcDim lonDim = dataFile.addDim("lon", NX);
//     NcDim latDim = dataFile.addDim("lat", NY);

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
//     countp.push_back(NX);
//     countp.push_back(NY);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         q_all.putVar(startp, countp, q[p]);
//     }
//     return;
// }

// void output_Qall(std::string dir,int n, double Q[6][NX][NY]) {
//     NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);       
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", NX);
//     NcDim yDim = dataFile.addDim("y", NY);
//     NcDim lonDim = dataFile.addDim("lon", NX);
//     NcDim latDim = dataFile.addDim("lat", NY);

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
//     countp.push_back(NX);
//     countp.push_back(NY);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         q_all.putVar(startp, countp, Q[p]);
//     }
//     return;
// }
