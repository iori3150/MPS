#include "settings.hpp"

#include <algorithm>
#include <fkYAML/node.hpp>
#include <format>
#include <fstream>
#include <iostream>

using std::format;

void Settings::load() {
    std::ifstream ifs("input/settings.yml");
    if (!ifs.is_open()) {
        std::cerr << format("Could not open setting file: {}", "input/settings.yml")
                  << std::endl;
        exit(-1);
    }

    fkyaml::node root = fkyaml::node::deserialize(ifs);

    try {
        dim              = root["dim"].get_value<int>();
        particleDistance = root["particleDistance"].get_value<double>();
        dt               = root["dt"].get_value<double>();
        finishTime       = root["finishTime"].get_value<double>();
        outputInterval   = root["outputInterval"].get_value<double>();
        cflCondition     = root["cflCondition"].get_value<double>();

        double xMin = root["domainRangeX"][0].get_value<double>();
        double xMax = root["domainRangeX"][1].get_value<double>();
        double yMin = root["domainRangeY"][0].get_value<double>();
        double yMax = root["domainRangeY"][1].get_value<double>();
        double zMin = root["domainRangeZ"][0].get_value<double>();
        double zMax = root["domainRangeZ"][1].get_value<double>();
        domain      = Domain(xMin, xMax, yMin, yMax, zMin, zMax);

        density            = root["density"].get_value<double>();
        kinematicViscosity = root["kinematicViscosity"].get_value<double>();

        gravity.x() = root["gravity"][0].get_value<double>();
        gravity.y() = root["gravity"][1].get_value<double>();
        gravity.z() = root["gravity"][2].get_value<double>();

        effectiveRadius.pressure =
            root["effectiveRadiusRatio"]["pressure"].get_value<double>() *
            particleDistance;
        effectiveRadius.viscosity =
            root["effectiveRadiusRatio"]["viscosity"].get_value<double>() *
            particleDistance;
        effectiveRadius.surfaceDetection =
            root["effectiveRadiusRatio"]["surfaceDetection"].get_value<double>() *
            particleDistance;
        effectiveRadius.max = std::max(
            {effectiveRadius.pressure,
             effectiveRadius.viscosity,
             effectiveRadius.surfaceDetection}
        );

        thresholdForSurfaceDetection =
            root["thresholdForSurfaceDetection"].get_value<double>();

        compressibility = root["compressibility"].get_value<double>();
        relaxationCoefficientForPressure =
            root["relaxationCoefficientForPressure"].get_value<double>();

        collisionDistance =
            root["collisionDistanceRatio"].get_value<double>() * particleDistance;
        coefficientOfRestitution = root["coefficientOfRestitution"].get_value<double>();

    } catch (const fkyaml::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }
}
