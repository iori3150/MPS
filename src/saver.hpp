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

private:
    std::filesystem::path outputDir;
};
