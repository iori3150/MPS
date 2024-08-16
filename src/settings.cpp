#include "settings.hpp"

#include <algorithm>
#include <fkYAML/node.hpp>
#include <format>
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>

using std::format;

void Settings::load(const std::filesystem::path& inputYamlPath) {
    std::ifstream ifs(inputYamlPath);
    if (!ifs.is_open()) {
        spdlog::error("Could not open setting file: " + inputYamlPath.string());
    }

    fkyaml::node root = fkyaml::node::deserialize(ifs);

    try {
        dim              = root["dim"].get_value<int>();
        particleDistance = root["particle distance"].get_value<double>();
        dt               = root["dt"].get_value<double>();
        finishTime       = root["finish time"].get_value<double>();
        outputInterval   = root["output interval"].get_value<double>();
        cflCondition     = root["CFL condition"].get_value<double>();

        double xMin = root["domain range X"][0].get_value<double>();
        double xMax = root["domain range X"][1].get_value<double>();
        double yMin = root["domain range Y"][0].get_value<double>();
        double yMax = root["domain range Y"][1].get_value<double>();
        double zMin = root["domain range Z"][0].get_value<double>();
        double zMax = root["domain range Z"][1].get_value<double>();
        domain      = Domain(xMin, xMax, yMin, yMax, zMin, zMax);

        density            = root["density"].get_value<double>();
        kinematicViscosity = root["kinematic viscosity"].get_value<double>();

        gravity.x() = root["gravity"][0].get_value<double>();
        gravity.y() = root["gravity"][1].get_value<double>();
        gravity.z() = root["gravity"][2].get_value<double>();

        effectiveRadius.pressure =
            root["effective radius ratio"]["pressure"].get_value<double>() *
            particleDistance;
        effectiveRadius.viscosity =
            root["effective radius ratio"]["viscosity"].get_value<double>() *
            particleDistance;
        effectiveRadius.surfaceDetection =
            root["effective radius ratio"]["surface detection"].get_value<double>() *
            particleDistance;
        effectiveRadius.max = std::max(
            {effectiveRadius.pressure,
             effectiveRadius.viscosity,
             effectiveRadius.surfaceDetection}
        );

        thresholdForSurfaceDetection =
            root["threshold for surface detection"].get_value<double>();

        compressibility = root["compressibility"].get_value<double>();
        relaxationCoefficientForPressure =
            root["relaxation coefficient for pressure"].get_value<double>();
        higherOrderSourceTerm.on =
            root["higher order source term"]["on"].get_value<bool>();
        if (higherOrderSourceTerm.on) {
            higherOrderSourceTerm.gamma =
                root["higher order source term"]["gamma"].get_value<double>();
        }

        collisionDistance =
            root["collision distance ratio"].get_value<double>() * particleDistance;
        coefficientOfRestitution = root["coefficient of restitution"].get_value<double>();

        inputCsvPath = root["inputCsvPath"].get_value<std::string>();
        inputCsvPath = inputYamlPath.parent_path() / inputCsvPath;

    } catch (const fkyaml::exception& e) {
        spdlog::error(e.what());
    }
}
