#include "simulation.hpp"

#include <argparse/argparse.hpp>
#include <cstdlib> // for std::exit
#include <filesystem>
#include <iostream>

using std::cout;
using std::endl;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("mps");

    program.add_argument("-s", "--setting")
        .required()
        .help("path to YAML file for settings.")
        .action([](const std::string& value) {
            if (!fs::exists(value)) {
                cout << "ERROR: Input YAML file cannot be found in specified path: " +
                            value
                     << endl;
                std::exit(EXIT_FAILURE);
            }
            cout << "Setting file: " << value << endl;
            return value;
        });

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        cout << err.what() << endl;
        std::exit(EXIT_FAILURE);
    }

    fs::path inputYamlPath = fs::path(program.get<std::string>("--setting"));

    Simulation simulation;
    simulation.run(inputYamlPath);

    return 0;
}
