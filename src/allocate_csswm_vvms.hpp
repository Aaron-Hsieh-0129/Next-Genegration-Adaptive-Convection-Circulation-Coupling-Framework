#include "../2DVVM/src/Declare.hpp"
#include "../CSSWM/src/construction.hpp"
#include "../CSSWM/src/define.hpp"

Config_VVM createConfig(const std::string& path, double addforcingtime, int CASE, double Kx, double Kz);
vvm**** allocate_and_initialize(int dim1, int dim2, int dim3);
Config_VVM**** allocate_and_initialize_config(int dim1, int dim2, int dim3);
void deallocate_config(Config_VVM**** array, int dim1, int dim2, int dim3);
void deallocate(vvm**** array, int dim1, int dim2, int dim3);