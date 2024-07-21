#pragma once

#include "particle.hpp"
#include "settings.hpp"

class MPS {
private:
    Settings settings;
    std::vector<Particle> particles;

public:
    MPS(Settings& settings, std::vector<Particle>& particles);

    void calGravity(std::vector<Particle>& particles);
};
