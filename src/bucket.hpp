#pragma once

#include "particle.hpp"
#include "settings.hpp"

#include <vector>

class Bucket {
private:
public:
    int num, numX, numY, numZ;
    double length;
    Domain domain;
    std::vector<int> next, first, last;

    Bucket() = default;
    Bucket(const double& reMax, const Domain& domain, const int& particleSize);
    void storeParticles(std::vector<Particle>& particles);
};
