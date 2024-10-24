#include <string>
#include <vector>
#include <map>

struct vvm_index {
    int p, i, j;
};

std::map<std::string, std::string> read_config(const std::string& filename);
std::vector<vvm_index> parse_int_tuples(const std::string& tuple_string);
bool is_value_in_vvm_index(const std::vector<vvm_index>& indices, int p, int i, int j);