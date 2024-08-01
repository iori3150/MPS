#include "simulation.hpp"

#include <filesystem>

int main() {
    std::filesystem::path inputYamlPath = "input/settings.yml";

    Simulation simulation;
    simulation.run(inputYamlPath);

    return 0;
}
