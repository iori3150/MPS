#include "particle.hpp"

#include <Eigen/Dense>

Particle::Particle(int id, ParticleType type, Eigen::Vector3d position, Eigen::Vector3d velocity) {
  this->id       = id;
  this->type     = type;
  this->position = position;
  this->velocity = velocity;
}