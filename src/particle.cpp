#include "particle.hpp"

#include <Eigen/Dense>

Particle::Particle(int id, ParticleType type, Eigen::Vector3d position) {
  this->id       = id;
  this->type     = type;
  this->position = position;
}