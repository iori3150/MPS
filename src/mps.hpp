#pragma once

#include "particle.hpp"
#include "settings.hpp"

#include <map>

class MPS {
private:
    struct {
        double gradient  = 0;
        double laplacian = 0;
    } n0;
    double lambda = 0;

    Settings settings;
    std::vector<Particle> particles;

    double weight(const double& dist, const double& re);

public:
    MPS(Settings& settings, std::vector<Particle>& particles);

    void calGravity(std::vector<Particle>& particles);
    void calViscosity(std::vector<Particle>& particles);
    void moveParticle(std::vector<Particle>& particles);
    void collision(std::vector<Particle>& particles);

    void calcNumberDensity(std::vector<Particle>& particles);
};
