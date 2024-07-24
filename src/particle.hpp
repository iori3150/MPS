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

enum class BoundaryCondition { Ignored, Surface, Inner };

class Particle {
private:
public:
    int id;
    ParticleType type;
    Eigen::Vector3d position, velocity;
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
    double pressure, minimumPressure = 0;
    double density;
    double numberDensityRatio = 0;

    std::vector<Neighbor> neighbors;

    BoundaryCondition boundaryCondition;
    bool isDirichletBoundaryConnected;

    Particle(
        const int& id,
        const ParticleType& type,
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const double& density
    );

    double inverseDensity() const;
};
