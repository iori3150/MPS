#pragma once

#include "particle.hpp"

#include <vector>

class Domain {
private:
public:
    struct {
        double min, max, length;
    } x, y, z;

    Domain() = default;
    Domain(
        const double& xMin,
        const double& xMax,
        const double& yMin,
        const double& yMax,
        const double& zMin,
        const double& zMax
    );
};

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
