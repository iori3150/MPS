#include "mps.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <iostream>

#define rep(i, a, b) for (int i = a; i < b; i++)
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0
#define DIRICHLET_BOUNDARY_IS_CONNECTED 1
#define DIRICHLET_BOUNDARY_IS_CHECKED 2

using std::cerr;
using std::cout;
using std::endl;

double weight(const double& dist, const double& re) {
    double w = 0.0;

    if (dist < re)
        w = (re / dist) - 1.0;

    return w;
}

MPS::MPS(const Settings& settings, std::vector<Particle>& particles) {
    this->settings  = settings;
    this->particles = particles;
    this->sourceTerm.resize(particles.size());
    this->coeffMatrix.resize(particles.size(), particles.size());
    this->flagForCheckingBoundaryCondition.resize(particles.size());
    this->bucket =
        Bucket(settings.effectiveRadius.max, settings.domain, particles.size());

    this->refValues.pressure = RefValues(
        settings.dim,
        settings.particleDistance,
        settings.effectiveRadius.pressure
    );
    this->refValues.viscosity = RefValues(
        settings.dim,
        settings.particleDistance,
        settings.effectiveRadius.viscosity
    );
    this->refValues.surfaceDetection = RefValues(
        settings.dim,
        settings.particleDistance,
        settings.effectiveRadius.surfaceDetection
    );
}

void MPS::stepForward(const bool isTimeToSave) {
    setNeighbors();
    calcGravity();
    calcViscosity();
    moveParticles();

    setNeighbors();
    collision();

    setNeighbors();
    calcPressure();
    calcPressureGradient();
    moveParticlesWithPressureGradient();

    if (isTimeToSave) {
        setNumberDensityForDisplay();
    }
}

RefValues::RefValues(
    const int& dim, const double& particleDistance, const double& effectiveRadius
) {
    int iZ_start, iZ_end;
    if (dim == 2) {
        iZ_start = 0;
        iZ_end   = 1;
    } else {
        iZ_start = -4;
        iZ_end   = 5;
    }
    for (int iX = -4; iX < 5; iX++) {
        for (int iY = -4; iY < 5; iY++) {
            for (int iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))
                    continue;

                Eigen::Vector3d r(
                    particleDistance * (double) iX,
                    particleDistance * (double) iY,
                    particleDistance * (double) iZ
                );
                double dist  = r.norm();
                double dist2 = dist * dist;

                initialNumberDensity += weight(dist, effectiveRadius);
                lambda += dist2 * weight(dist, effectiveRadius);
            }
        }
    }
    lambda /= initialNumberDensity;
}

void MPS::calcGravity() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.acceleration = settings.gravity;

        } else {
            pi.acceleration.setZero();
        }
    }
}

void MPS::calcViscosity() {
    const double& re     = settings.effectiveRadius.viscosity;
    const double& n0     = refValues.viscosity.initialNumberDensity;
    const double& lambda = refValues.viscosity.lambda;

    double a = (settings.kinematicViscosity) * (2.0 * settings.dim) / (n0 * lambda);

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type != ParticleType::Fluid)
            continue;

        Eigen::Vector3d viscosity_term = Eigen::Vector3d::Zero();

        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double dist  = neighbor.distance;

            if (dist < re) {
                viscosity_term += (pj.velocity - pi.velocity) * weight(dist, re);
            }
        }

        viscosity_term *= a;
        pi.acceleration += viscosity_term;
    }
}

void MPS::moveParticles() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.velocity += pi.acceleration * settings.dt;
            pi.position += pi.velocity * settings.dt;
        }
        pi.acceleration.setZero();
    }
}

void MPS::collision() {
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

void MPS::calcPressure() {
    setBoundaryCondition();
    setSourceTerm();
    setMatrix();
    solvePoissonEquation();
    removeNegativePressure();
    setMinimumPressure();
}

void MPS::setBoundaryCondition() {
    double n0   = refValues.surfaceDetection.initialNumberDensity;
    double beta = settings.thresholdForSurfaceDetection;

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost || pi.type == ParticleType::DummyWall) {
            pi.boundaryCondition = BoundaryCondition::GhostOrDummy;

        } else if (getNumberDensity(pi, settings.effectiveRadius.surfaceDetection) < beta * n0) {
            pi.boundaryCondition = BoundaryCondition::Surface;

        } else {
            pi.boundaryCondition = BoundaryCondition::Inner;
        }
    }
}

void MPS::setSourceTerm() {
    const double re    = settings.effectiveRadius.pressure;
    const double n0    = refValues.pressure.initialNumberDensity;
    const double gamma = settings.relaxationCoefficientForPressure;

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::Inner) {
            sourceTerm[pi.id] = gamma * (1.0 / (settings.dt * settings.dt)) *
                                ((getNumberDensity(pi, re) - n0) / n0);

        } else {
            sourceTerm[pi.id] = 0.0;
        }
    }
}

void MPS::setMatrix() {
    coeffMatrix.setZero();

    const double re     = settings.effectiveRadius.pressure;
    const double n0     = refValues.pressure.initialNumberDensity;
    const double lambda = refValues.pressure.lambda;
    const double a      = 2.0 * settings.dim / (n0 * lambda);
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition != BoundaryCondition::Inner)
            continue;

        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re) {
                double coef_ij = a * weight(dist, re) / pi.density;

                coeffMatrix(pi.id, pj.id) = (-1.0 * coef_ij);
                coeffMatrix(pi.id, pi.id) += coef_ij;
            }
        }

        coeffMatrix(pi.id, pi.id) +=
            (settings.compressibility) / (settings.dt * settings.dt);
    }

    exceptionalProcessingForBoundaryCondition();
}

void MPS::exceptionalProcessingForBoundaryCondition() {
    // If tere is no Dirichlet boundary condition on the fluid,
    // increase the diagonal terms of the matrix for an exception.
    // This allows us to solve the matrix without Dirichlet boundary conditions.
    checkBoundaryCondition();
    increaseDiagonalTerm();
}

void MPS::solvePoissonEquation() {
    Eigen::SparseMatrix<double> A = coeffMatrix.sparseView();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd pressures = solver.solve(sourceTerm);
    for (auto& pi : particles) {
        pi.pressure = pressures[pi.id];
    }
}

void MPS::checkBoundaryCondition() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::GhostOrDummy) {
            flagForCheckingBoundaryCondition[pi.id] = -1;

        } else if (pi.boundaryCondition == BoundaryCondition::Surface) {
            flagForCheckingBoundaryCondition[pi.id] = DIRICHLET_BOUNDARY_IS_CONNECTED;

        } else {
            flagForCheckingBoundaryCondition[pi.id] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
        }
    }

    int count;
    while (true) {
        count = 0;

        rep(i, 0, particles.size()) {
            if (flagForCheckingBoundaryCondition[i] != DIRICHLET_BOUNDARY_IS_CONNECTED)
                continue;

            for (auto& neighbor : particles[i].neighbors) {
                const Particle& pj = particles[neighbor.id];
                const double& dist = neighbor.distance;
                const double& re   = settings.effectiveRadius.pressure;

                if (flagForCheckingBoundaryCondition[pj.id] !=
                    DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                    continue;

                if (dist < re)
                    flagForCheckingBoundaryCondition[pj.id] =
                        DIRICHLET_BOUNDARY_IS_CONNECTED;
            }

            flagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
            count++;
        }

        if (count == 0)
            break;
    }

#pragma omp parallel for
    rep(i, 0, particles.size()) {
        if (flagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
#pragma omp critical
            {
                cerr << "WARNING: There is no dirichlet boundary condition for particle "
                     << i << endl;
            }
        }
    }
}

void MPS::increaseDiagonalTerm() {
#pragma omp parallel for
    rep(i, 0, particles.size()) {
        if (flagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
            coeffMatrix(i, i) *= 2.0;
        }
    }
}

void MPS::removeNegativePressure() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.pressure < 0.0) {
            pi.pressure = 0.0;
        }
    }
}

void MPS::setMinimumPressure() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::GhostOrDummy)
            continue;

        pi.minimumPressure = pi.pressure;

        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;
            const double& re   = settings.effectiveRadius.pressure;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re) {
                if (pi.minimumPressure > pj.pressure)
                    pi.minimumPressure = pj.pressure;
            }
        }
    }
}

void MPS::calcPressureGradient() {
    const double& re = settings.effectiveRadius.pressure;
    const double& n0 = refValues.pressure.initialNumberDensity;

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type != ParticleType::Fluid)
            continue;

        Eigen::Vector3d gradient = Eigen::Vector3d::Zero();

        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re) {
                gradient += (pj.position - pi.position) *
                            (pj.pressure - pi.minimumPressure) * weight(dist, re) /
                            (dist * dist);
            }
        }

        gradient *= settings.dim / n0;
        pi.acceleration = -1.0 * gradient / pi.density;
    }
}

void MPS::moveParticlesWithPressureGradient() {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.velocity += pi.acceleration * settings.dt;
            pi.position += pi.acceleration * settings.dt * settings.dt;
        }

        pi.acceleration.Zero();
    }
}

double MPS::getCourantNumber() {
    double maxCourantNumber = 0.0;

    for (auto& pi : particles) {
        if (pi.type != ParticleType::Fluid)
            continue;

        double courantNumeber =
            (pi.velocity.norm() * settings.dt) / settings.particleDistance;
        maxCourantNumber = std::max(maxCourantNumber, courantNumeber);
    }

    if (maxCourantNumber > settings.cflCondition) {
        cerr << "ERROR: Courant number is larger than CFL condition. Courant = "
             << maxCourantNumber << endl;
    }

    return maxCourantNumber;
}

double MPS::getNumberDensity(const Particle& pi, const double& re) {
    double numberDensity = 0.0;
    for (auto& neighbor : pi.neighbors) {
        const double& dist = neighbor.distance;

        if (dist < re) {
            numberDensity += weight(dist, re);
        }
    }
    return numberDensity;
}

void MPS::setNeighbors() {
    bucket.storeParticles(particles);

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost)
            continue;

        pi.neighbors.clear();

        int ix = int((pi.position.x() - bucket.domain.x.min) / bucket.length) + 1;
        int iy = int((pi.position.y() - bucket.domain.y.min) / bucket.length) + 1;
        int iz = int((pi.position.z() - bucket.domain.z.min) / bucket.length) + 1;

        for (int jx = ix - 1; jx <= ix + 1; jx++) {
            for (int jy = iy - 1; jy <= iy + 1; jy++) {
                for (int jz = iz - 1; jz <= iz + 1; jz++) {
                    int jBucket = jx + jy * bucket.numX + jz * bucket.numX * bucket.numY;
                    int j       = bucket.first[jBucket];

                    while (j != -1) {
                        Particle& pj = particles[j];

                        double dist = (pj.position - pi.position).norm();
                        if (j != pi.id && dist < settings.effectiveRadius.max) {
                            pi.neighbors.emplace_back(j, dist);
                        }

                        j = bucket.next[j];
                    }
                }
            }
        }
    }
}

// Set the number density at the final position so that users can visually comprehend how
// dense/sparse the particles are.
// This function should only be called when it's time to save, because setting neighbors
// and number density is actually unnecessary for the simulation.
void MPS::setNumberDensityForDisplay() {
    setNeighbors();

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::Surface) {
            pi.numberDensityRatio = 0.0;

        } else {
            pi.numberDensityRatio =
                getNumberDensity(pi, settings.effectiveRadius.pressure) /
                refValues.pressure.initialNumberDensity;
        }
    }
}
