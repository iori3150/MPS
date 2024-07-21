#include "mps.hpp"

#define rep(i, a, b) for (int i = a; i < b; i++)

MPS::MPS(Settings& settings, std::vector<Particle>& particles) {
    this->settings  = settings;
    this->particles = particles;

    int iZ_start = -4;
    int iZ_end   = 5;
    if (settings.dim == 2) {
        iZ_start = 0;
        iZ_end   = 1;
    }
    for (int iX = -4; iX < 5; iX++) {
        for (int iY = -4; iY < 5; iY++) {
            for (int iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))
                    continue;

                Eigen::Vector3d r(
                    settings.particleDistance * (double) iX,
                    settings.particleDistance * (double) iY,
                    settings.particleDistance * (double) iZ
                );
                double dist  = r.norm();
                double dist2 = dist * dist;

                n0.gradient += weight(dist, settings.re.gradient);
                n0.laplacian += weight(dist, settings.re.laplacian);
                lambda += dist2 * weight(dist, settings.re.laplacian);
            }
        }
    }
    lambda /= n0.laplacian;
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

void MPS::calViscosity(std::vector<Particle>& particles) {
    const double& n0 = this->n0.laplacian;
    const double& re = settings.re.laplacian;

    double A = (settings.kinematicViscosity) * (2.0 * settings.dim) / (n0 * lambda);

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type != ParticleType::Fluid)
            continue;

        Eigen::Vector3d viscosity_term = Eigen::Vector3d::Zero();

        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double dist  = neighbor.distance;

            if (dist < settings.re.laplacian) {
                viscosity_term += (pj.velocity - pi.velocity) * weight(dist, re);
            }
        }

        viscosity_term *= A;
        pi.acceleration += viscosity_term;
    }
}

double MPS::weight(const double& dist, const double& re) {
    double w = 0.0;

    if (dist < re)
        w = (re / dist) - 1.0;

    return w;
}
