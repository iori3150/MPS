#include "simulation.hpp"

#include <cstdlib> // for std::exit
#include <filesystem>
#include <iostream>

int main() {
    std::filesystem::path inputYamlPath = "input/settings.yml";
    if (!std::filesystem::exists(inputYamlPath)) {
        std::cout << "Input YAML file cannot be found in specified path: " +
                         inputYamlPath.string()
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Simulation simulation;
    simulation.run(inputYamlPath);

    return 0;
}
