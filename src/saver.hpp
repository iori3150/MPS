#pragma once

#include "particle.hpp"

#include <filesystem>
#include <vector>

class Saver {
public:
    Saver() = default;
    Saver(const std::filesystem::path& outputDir);

    void
    save(const std::vector<Particle>& particles, const double& time, const int& fileNum);

private:
    std::filesystem::path outputDir;
};
