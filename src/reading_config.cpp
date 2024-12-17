#include "reading_config.hpp"
#include <fstream>
#include <sstream>


std::map<std::string, std::string> read_config(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, std::string> config;
    std::string line;

    // Read each line from the file
    while (std::getline(file, line)) {
        // Ignore everything after a '#' character
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Remove leading and trailing whitespace (optional)
        line.erase(0, line.find_first_not_of(" \t")); // Leading
        line.erase(line.find_last_not_of(" \t") + 1); // Trailing

        if (line.empty()) {
            continue; // Skip empty lines
        }

        std::istringstream is_line(line);
        std::string key;

        // Extract the key before '='
        if (std::getline(is_line, key, '=')) {
            std::string value;

            // Remove leading/trailing whitespace from the key
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);

            // Extract the value after '='
            if (std::getline(is_line, value)) {
                // Remove leading/trailing whitespace from the value
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                // Store the key-value pair in the map
                config[key] = value;
            }
        }
    }

    return config; // Return the map with all key-value pairs
}

// Function to parse a list of integer tuples in the format [(p, i, j), (p2, i2, j2), ...]
std::vector<vvm_index> parse_int_tuples(const std::string& tuple_string) {
    std::vector<vvm_index> indices;
    int p = 0, i = 0, j = 0;
    bool in_number = false;
    int num = 0;
    int count = 0;

    // Loop through the string to manually extract numbers
    for (size_t idx = 0; idx < tuple_string.size(); ++idx) {
        char c = tuple_string[idx];

        // If character is a digit, update the current number
        if (isdigit(c)) {
            num = num * 10 + (c - '0');
            in_number = true;
        }
        // If character is not a digit but we just finished reading a number
        else if (in_number) {
            if (count == 0) p = num;
            else if (count == 1) i = num;
            else if (count == 2) j = num;

            // Reset for the next number
            num = 0;
            in_number = false;
            count++;
        }

        // If we have three numbers (a full tuple), add to vector and reset for the next one
        if (count == 3) {
            indices.push_back(vvm_index{p, i, j}); // Add the struct to the array
            count = 0; // Reset the counter for the next tuple
        }
    }

    return indices;
}

bool is_value_in_vvm_index(const std::vector<vvm_index>& indices, int p, int i, int j) {
    for (const auto& idx : indices) {
        if (idx.p == p && idx.i == i && idx.j == j) {
            return true;
        }
    }
    return false;
}