#pragma once

#include "particle.hpp"

#include <filesystem>
#include <vector>

class Saver {
public:
    int numberOfFiles = 0;

    Saver() = default;
    Saver(const std::filesystem::path& outputDir);

    void save(const std::vector<Particle>& particles, const double& time);
    void toCsv(const std::vector<Particle>& particles, const double& time);
    void toVtu(const std::vector<Particle>& particles, const double& time);
    void dataArrayBegin(
        std::ofstream& outFile,
        const std::string& numberOfComponents,
        const std::string& type,
        const std::string& name
    );
    void dataArrayEnd(std::ofstream& outFile);

private:
    std::filesystem::path outputDir;
};
