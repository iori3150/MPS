#pragma once

#include "particle.hpp"
#include "settings.hpp"

#include <map>

class MPS {
private:
    std::map<std::string, double> n0 = {{"gradient", 0}, {"laplacian", 0}};
    double lambda                    = 0;

    Settings settings;
    std::vector<Particle> particles;

    double weight(const double& dist, const double& re);

public:
    MPS(Settings& settings, std::vector<Particle>& particles);

    void calGravity(std::vector<Particle>& particles);
    void calViscosity(std::vector<Particle>& particles);
};
