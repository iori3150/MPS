#pragma once

#include "Eigen/Dense"

enum class ParticleType { Ghost = -1, Fluid = 0, Wall = 1, DummyWall = 2 };

class Neighbor {
private:
public:
    int id;
    double distance;

    Neighbor(int id, double distance) {
        this->id       = id;
        this->distance = distance;
    };
};

class Particle {
private:
public:
    int id;
    ParticleType type;
    Eigen::Vector3d position, velocity;
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
    double pressure;
    double numberDensity = 0;
    double density;

    std::vector<Neighbor> neighbors;

    Particle(
        const int& id,
        const ParticleType& type,
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const double& pressure,
        const double& density
    );

    double inverseDensity() const;
};
