#include "simulation.hpp"

#include <filesystem>
#include <iostream>

int main() {
    std::filesystem::path inputYamlPath = "input/settings.yml";
    if (!std::filesystem::exists(inputYamlPath)) {
        std::cout << "Input YAML file cannot be found in specified path: "
                  << inputYamlPath << std::endl;
        exit(-1);
    }

    Simulation simulation;
    simulation.run(inputYamlPath);

    return 0;
}
