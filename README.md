# Next Genegration Adaptive Convection-Circulation Coupling Framework


This is a framework coupled a cubed-sphere shallow water model and a 2D vector vorticity cloud resolving model

It's a developing project so the interface is not constructed completely. If you want to modify the logic inside the coupler or the two models, please check the Guidance for Code Modification.

- [Next Genegration Adaptive Convection-Circulation Coupling Framework](#next-genegration-adaptive-convection-circulation-coupling-framework)
  - [Prerequisite](#prerequisite)
  - [How to Use](#how-to-use)
    - [Guidance for Code Modification](#guidance-for-code-modification)

## Prerequisite

- C++ compiler (higher than C++11)
- CMake (higher than 3.0.0) (You can create your own Makefile by translating the CMakefile.txt if you don't want to use CMake)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx) [optional]
- PETSc [optional]
- Eigen (this has already be installed in include folder) [optional]

  
## How to Use

1. Clone the project using

   ```
   git clone --recurse-submodules https://github.com/Aaron-Hsieh-0129/Next-Genegration-Adaptive-Convection-Circulation-Coupling-Framework.git
   ```

2. You are able to run the model by running the command under the project folder

    ```
    sh TMIF.sh
    ```

or you can use your own command by referencing the command in `TMIF.sh`


### Guidance for Code Modification
- This projects are composed of a coupler in `src/main.cpp` and two modules `CSSWM` and `2DVVM`. 
- If you want to tune the input variables of `CSSWM` and `2DVVM`, it's recommended to change the configurations of them in `src/main.cpp` rather than modify the settings in `CSSWM` and `2DVVM`
- There are three variables that can be changed in `config.txt` , and more options will be appended in the future to make the user no need to modify `src/main.cpp`
    ```
    OUTPUTPATH=/work/aaron900129/NGC3F/DATA/test/
    SEED=0
    COUPLETIME=600
    Bubble_p_i_j=[(1,46,47), (1,47,47), (1,48,47)]
    NotBubble_p_i_j=[(1,44,47), (1,45,47), (1,49,47), (1,50,47)]
    ```
