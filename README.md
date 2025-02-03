# Next Genegration Adaptive Convection-Circulation Coupling Framework


This is a framework coupled a cubed-sphere shallow water model and a 2D vector vorticity cloud resolving model

It's a developing project so the interface is not constructed completely. If you want to modify the logic inside the coupler or the two models, please check the Guidance for Code Modification.

- [Next Genegration Adaptive Convection-Circulation Coupling Framework](#next-genegration-adaptive-convection-circulation-coupling-framework)
  - [Prerequisite](#prerequisite)
  - [How to Use](#how-to-use)
    - [Guidance for Code Modification](#guidance-for-code-modification)
    - [Results Quickview](#results-quickview)

## Prerequisite

- C++ compiler (higher than C++11)
- CMake (higher than 3.0.0) (You can create your own Makefile by translating the CMakefile.txt if you don't want to use CMake)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx) [optional]
- PETSc [optional]
- Eigen (this has already be installed in include folder) [optional]

  
## How to Use

1. Clone the project using

   ```
   git clone --recurse-submodules https://github.com/Aaron-Hsieh-0129/Next-Genegration-Adaptive-Convection-Circulation-Coupling-Framework.git NextACC
   ```

2. You are able to run the model by running the command under the project folder

    ```
    sh TMIF.sh
    ```

or you can use your own command by referencing the command in `TMIF.sh`

Noted that you might encounter some errors related to the compiler or the libraries, please modify the paths of the compiler and the libraries path in `./CMakeLists.txt` to the place you install them. 


### Guidance for Code Modification
- This projects are composed of a coupler in `src/main.cpp` and two modules `CSSWM` and `2DVVM`. 
- If you want to tune the input variables of `CSSWM` and `2DVVM`, it's recommended to change the configurations of them in `src/main.cpp` rather than modify the settings in `CSSWM` and `2DVVM`
- There are some variables that can be changed in `./config.txt` to specify in CSSWM and 2DVVM, please try to modify here first. If you want some more advanced functions, then go into the 2DVVM and CSSWM folder to modify them. 
    ```
    OUTPUTPATH=~/NextACC/DATA/test
    SEED=0                                          # random seed for initial perturbation in the CRM bottom
    COUPLETIME=600                                  # Coupling time for NextGCC [s]
    Bubble_p_i_j=[]    # CRMs with bubble inside
    NotBubble_p_i_j=[(0,42,50),(0,44,50),(0,46,50),(0,49,47),(0,51,47)]             # CRMs with nothing inside
    BubbleCase=1                                    # Case0: Nothing, Case1: Bubble, Case2: Bubble+wind shear
    CSSWM_GRAVITY=0.2391                               # gravity wave speed for CSSWM

    CSSWM_DT=200
    CSSWM_TIMEEND=10000                             # Integration Time [s]
    CSSWM_OUTPUTSTEP=1                              # Output frequency
    CSSWM_DIFFUSION_KX=200000
    CSSWM_DIFFUSION_KY=200000
    CSSWM_DIFFUSION_TS=0.06
    CSSWM_ADDFORCING_TIME=6000                     # If the user specifies adding forcing, the adding time can be specified here
    CSSWM_H_NUDGE_TIME=0                           # CSSWM h nudging time scale, if it is 0, the nudge will be closed.


    VVM_XRANGE=20000                               # Domain for x [m]
    VVM_ZRANGE=20000                                # Domain for z [m]
    VVM_DT=3
    VVM_DX=200
    VVM_DZ=200                                      # Should be same as dx
    VVM_TIMEEND=10000                               # Integration Time [s]
    VVM_OUTPUTSTEP=50                               # Output frequency
    ```


### Results Quickview
- If you have GrADs installed, you can go to the data output path for CSSWM or 2DVVMs. There will be .ctl file created inside and you can directly use them to quickly view the results.
- Python visualization is also available and the visualize code will be put under `./scripts` folder.

