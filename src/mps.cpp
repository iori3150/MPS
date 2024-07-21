#include "mps.hpp"

MPS::MPS(Settings& settings, std::vector<Particle>& particles) {
    this->settings  = settings;
    this->particles = particles;
}

void MPS::calGravity(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.acceleration = settings.gravity;

        } else {
            pi.acceleration.setZero();
        }
    }
}
