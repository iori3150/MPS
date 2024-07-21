#include "particle.hpp"

#include <Eigen/Dense>

Particle::Particle(
    const int& id,
    const ParticleType& type,
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    const double& pressure,
    const double& density
) {
    this->id       = id;
    this->type     = type;
    this->position = position;
    this->velocity = velocity;
    this->pressure = pressure;
    this->density  = density;
}

double Particle::inverseDensity() const {
    switch (this->type) {
    case ParticleType::Fluid:
        return 1 / this->density;

    case ParticleType::Wall:
    case ParticleType::DummyWall:
    case ParticleType::Ghost:
        return 0;
    }
}
