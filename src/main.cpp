#include "simulation.hpp"
#include "utilities.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>

int main() {
    std::filesystem::path inputYamlPath = "input/settings.yml";
    if (!std::filesystem::exists(inputYamlPath)) {
        exitWithError(
            "Input YAML file cannot be found in specified path: " + inputYamlPath.string()
        );
    }

    Simulation simulation;
    simulation.run(inputYamlPath);

    return 0;
}
