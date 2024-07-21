#include "mps.hpp"

#define rep(i, a, b) for (int i = a; i < b; i++)

MPS::MPS(const Settings& settings, const int& numberOfParticles) {
    this->settings = settings;
    this->sourceTerm.resize(numberOfParticles);

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

                n0.numberDensity += weight(dist, settings.re.numberDensity);
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

void MPS::moveParticle(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.velocity += pi.acceleration * settings.dt;
            pi.position += pi.velocity * settings.dt;
        }
        pi.acceleration.setZero();
    }
}

void MPS::collision(std::vector<Particle>& particles) {
    for (auto& pi : particles) {
        if (pi.type != ParticleType::Fluid)
            continue;

        for (auto& neighbor : pi.neighbors) {
            Particle& pj = particles[neighbor.id];
            double& dist = neighbor.distance;

            if (pj.type == ParticleType::Fluid && pj.id >= pi.id)
                continue;

            if (dist < settings.collisionDistance) {
                double invMassi    = pi.inverseDensity();
                double invMassj    = pj.inverseDensity();
                double reducedMass = 1.0 / (invMassi + invMassj);

                Eigen::Vector3d normal  = (pj.position - pi.position).normalized();
                double relativeVelocity = (pj.velocity - pi.velocity).dot(normal);
                if (relativeVelocity < 0.0) {
                    double impulse = -(1.0 + settings.coefficientOfRestitution) *
                                     relativeVelocity * reducedMass;
                    pi.velocity -= impulse * invMassi * normal;
                    pj.velocity += impulse * invMassj * normal;

                    double depth           = settings.collisionDistance - dist;
                    double positionImpulse = depth * reducedMass;
                    pi.position -= positionImpulse * invMassi * normal;
                    pj.position += positionImpulse * invMassj * normal;

                    // cerr << "WARNING: Collision between particles " << pi.id << " and "
                    // << pj.id << " occurred."
                    // << endl;
                }
            }
        }
    }
}

void MPS::calcNumberDensity(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost)
            continue;

        pi.numberDensity = 0.0;

        for (auto& neighbor : pi.neighbors) {
            const double& dist = neighbor.distance;
            const double& re   = settings.re.numberDensity;

            if (dist < re) {
                pi.numberDensity += weight(dist, re);
            }
        }
    }
}

void MPS::setBoundaryCondition(std::vector<Particle>& particles) {
    double n0   = this->n0.numberDensity;
    double beta = settings.thresholdForSurfaceDetection;

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost || pi.type == ParticleType::DummyWall) {
            pi.boundaryCondition = BoundaryCondition::GhostOrDummy;

        } else if (pi.numberDensity < beta * n0) {
            pi.boundaryCondition = BoundaryCondition::Surface;

        } else {
            pi.boundaryCondition = BoundaryCondition::Inner;
        }
    }
}

double MPS::weight(const double& dist, const double& re) {
    double w = 0.0;

    if (dist < re)
        w = (re / dist) - 1.0;

    return w;
}
