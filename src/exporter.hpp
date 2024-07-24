#pragma once

#include "particle.hpp"

#include <filesystem>
#include <vector>

class Exporter {
public:
    void toCsv(
        const std::filesystem::path& outputFilePath,
        const std::vector<Particle>& particles,
        const double& time
    );
    void toVtu(
        const std::filesystem::path& outputFilePath,
        const std::vector<Particle>& particles,
        const double& time
    );

private:
    void dataArrayBegin(
        std::ofstream& outFile,
        const std::string& numberOfComponents,
        const std::string& type,
        const std::string& name
    );
    void dataArrayEnd(std::ofstream& outFile);
};
