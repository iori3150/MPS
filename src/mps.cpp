#include "mps.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <csv.hpp>
#include <iostream>
#include <queue>
#include <spdlog/spdlog.h>

double weight(const double& dist, const double& re) {
    double w = 0.0;

    if (dist < re)
        w = (re / dist) - 1.0;

    return w;
}

MPS::MPS(const std::filesystem::path& inputYamlPath) {
    settings.load(inputYamlPath);

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

double MPS::importInitialCondition() {
    if (!std::filesystem::exists(settings.inputCsvPath)) {
        spdlog::error(
            "Input file does not exist in the specified path: {}",
            settings.inputCsvPath.string()
        );
    }

    int particleDataHeaderRow = 3;
    double initialTime;

    // Set up CSV format for meta data
    csv::CSVFormat metaDataFormat;
    metaDataFormat.no_header();
    csv::CSVReader metaDataReader(settings.inputCsvPath.string(), metaDataFormat);

    // Read meta data
    int rowNumber = 1;
    for (auto& row : metaDataReader) {
        // Get the time from the first row
        if (rowNumber == 1)
            initialTime = row[1].get<double>();

        // Stop reading after the row before the particle data header
        if (rowNumber == particleDataHeaderRow - 1)
            break;

        rowNumber++;
    }

    // Set up CSV format for particle data
    csv::CSVFormat particleDataFormat;
    particleDataFormat.header_row(particleDataHeaderRow - 1);
    csv::CSVReader particleDataReader(settings.inputCsvPath.string(), particleDataFormat);

    // Read particle data and create Particle objects
    for (auto& row : particleDataReader) {
        int id      = row["ID"].get<int>();
        double type = row["Type"].get<int>();
        double x    = row["Position.x (m)"].get<double>();
        double y    = row["Position.y (m)"].get<double>();
        double z    = row["Position.z (m)"].get<double>();
        double u    = row["Velocity.x (m/s)"].get<double>();
        double v    = row["Velocity.y (m/s)"].get<double>();
        double w    = row["Velocity.z (m/s)"].get<double>();

        particles.push_back(Particle(
            id,
            static_cast<ParticleType>(type),
            Eigen::Vector3d(x, y, z),
            Eigen::Vector3d(u, v, w),
            settings.density
        ));
    }

    return initialTime;
}

double MPS::loadInitialState() {
    double initialTime = importInitialCondition();

    this->sourceTerm.resize(particles.size());
    this->coefficientMatrix.resize(particles.size(), particles.size());
    this->bucket =
        Bucket(settings.effectiveRadius.max, settings.domain, particles.size());

    // Set the particle number density at the beginning so that we can check if the
    // initial particle placement is appropriate.
    setNumberDensity();

    return initialTime;
}

double MPS::stepForward(const bool isTimeToExport) {
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

    if (isTimeToExport || settings.higherOrderSourceTerm.on) {
        setNumberDensity();
    }

    return getCourantNumber();
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
            checkBoundaryViolation(pi);
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
                }
            }
        }
    }
}

void MPS::calcPressure() {
    setBoundaryCondition();
    setSourceTerm();
    setMatrix();
    ensureDirichletBoundaryConnection();
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
            pi.boundaryCondition = BoundaryCondition::Ignored;

        } else if (getNumberDensity(pi, settings.effectiveRadius.surfaceDetection) < beta * n0) {
            pi.boundaryCondition = BoundaryCondition::Surface;

        } else {
            pi.boundaryCondition = BoundaryCondition::Inner;
        }
    }
}

void MPS::setSourceTerm() {
    const double re = settings.effectiveRadius.pressure;
    const double n0 = refValues.pressure.initialNumberDensity;
    const double dt = settings.dt;

#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::Inner) {
            if (settings.higherOrderSourceTerm.on) {
                const double gamma = settings.higherOrderSourceTerm.gamma;
                sourceTerm[pi.id] =
                    gamma * (1.0 / (dt * dt)) * (pi.numberDensityRatio - 1.0);

                for (auto& neighbor : pi.neighbors) {
                    const Particle& pj = particles[neighbor.id];

                    if (pj.boundaryCondition == BoundaryCondition::Ignored)
                        continue;

                    const Eigen::Vector3d& r_ij = pj.position - pi.position;
                    const Eigen::Vector3d& u_ij = pj.velocity - pi.velocity;
                    const double& dist          = neighbor.distance;

                    if (dist < re) {
                        sourceTerm[pi.id] -= (1.0 / (n0 * dt)) * re * r_ij.dot(u_ij) /
                                             (dist * dist * dist);
                    }
                }

            } else {
                const double gamma = settings.relaxationCoefficientForPressure;
                const double n     = getNumberDensity(pi, re);

                sourceTerm[pi.id] = gamma * (1.0 / (dt * dt)) * (n - n0) / n0;
            }

        } else {
            sourceTerm[pi.id] = 0.0;
        }
    }
}

void MPS::setMatrix() {
    std::vector<Eigen::Triplet<double>> matrixTriplets;

    const double re     = settings.effectiveRadius.pressure;
    const double n0     = refValues.pressure.initialNumberDensity;
    const double lambda = refValues.pressure.lambda;
    const double D      = settings.dim;

    for (auto& pi : particles) {
        if (pi.boundaryCondition != BoundaryCondition::Inner)
            continue;

        double coefficient_ii = 0.0;
        for (auto& neighbor : pi.neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.boundaryCondition == BoundaryCondition::Ignored)
                continue;

            if (dist < re) {
                double a = (2.0 * D / (lambda * n0)) * weight(dist, re) / pi.density;

                coefficient_ii += a;
                if (pj.boundaryCondition == BoundaryCondition::Inner) {
                    matrixTriplets.emplace_back(pi.id, pj.id, -a);
                }
            }
        }

        if (settings.quasiCompressibility.on) {
            const double alpha = settings.quasiCompressibility.compressibility;
            coefficient_ii += alpha / (settings.dt * settings.dt);
        }
        matrixTriplets.emplace_back(pi.id, pi.id, coefficient_ii);
    }

    coefficientMatrix.setFromTriplets(matrixTriplets.begin(), matrixTriplets.end());
}

void MPS::ensureDirichletBoundaryConnection() {
    for (auto& pi : particles) {
        if (pi.boundaryCondition == BoundaryCondition::Surface) {
            pi.isDirichletBoundaryConnected = true;
        } else {
            pi.isDirichletBoundaryConnected = false;
        }
    }

    // Beginning from Surface particles, perform BFS (Breath First Search) to check
    // Dirichlet boundary connection
    for (auto& p : particles) {
        if (p.boundaryCondition != BoundaryCondition::Surface)
            continue;

        // Create a queue for BFS and push the id of the Surface particle
        std::queue<int> queue;
        queue.push(p.id);

        // Loop through the queue until it's empty
        while (!queue.empty()) {
            const Particle& pi = particles[queue.front()];
            queue.pop();

            // Look through all neighbors of particle i
            for (auto& neighbor : pi.neighbors) {
                Particle& pj       = particles[neighbor.id];
                const double& dist = neighbor.distance;
                const double& re   = settings.effectiveRadius.pressure;

                // If the neighbor is witin the effective radius and not yet marked, mark
                // it as connected and add to the queue
                if (dist < re && !pj.isDirichletBoundaryConnected) {
                    pj.isDirichletBoundaryConnected = true;
                    queue.push(pj.id);
                }
            }
        }
    }

    for (auto& pi : particles) {
        if (!pi.isDirichletBoundaryConnected &&
            pi.boundaryCondition == BoundaryCondition::Inner) {
            spdlog::debug(
                "There is no Dirichlet boundary condition connected to the particle "
                "(id={}).",
                std::to_string(pi.id)
            );
            debugLogCount++;

            coefficientMatrix.coeffRef(pi.id, pi.id) *= 2.0;
        }
    }
}

void MPS::solvePoissonEquation() {
    if (coefficientMatrix.isApprox(coefficientMatrix.transpose())) {
        // std::cout << "Matrix is symmetric." << std::endl;
    } else {
        for (int i = 0; i < coefficientMatrix.rows(); ++i) {
            for (int j = i + 1; j < coefficientMatrix.cols(); ++j) {
                if (std::abs(
                        coefficientMatrix.coeff(i, j) - coefficientMatrix.coeff(j, i)
                    ) > 1.0e-10) {
                    std::cout << "Matrix is not symmetric at (" << i << ", " << j << "): "
                              << "A(" << i << ", " << j
                              << ") = " << coefficientMatrix.coeff(i, j) << ", A(" << j
                              << ", " << i << ") = " << coefficientMatrix.coeff(j, i)
                              << std::endl;
                    std::cout << i << ": type=" << static_cast<int>(particles[i].type)
                              << " BC="
                              << static_cast<int>(particles[i].boundaryCondition)
                              << std::endl;
                    std::cout << j << ": type=" << static_cast<int>(particles[j].type)
                              << " BC="
                              << static_cast<int>(particles[j].boundaryCondition)
                              << std::endl;
                    std::exit(1);
                }
            }
        }
    }

    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.analyzePattern(coefficientMatrix); // パターン解析
    // solver.factorize(coefficientMatrix);      // LU分解
    // if (solver.info() != Eigen::Success) {
    //     spdlog::error("Factorization failed.");
    // }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(coefficientMatrix);

    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    // solver.compute(coefficientMatrix);
    Eigen::VectorXd pressures = solver.solve(sourceTerm);
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error() << std::endl;
    pressures = solver.solve(sourceTerm);

    if (solver.info() != Eigen::Success) {
        // spdlog::error("Pressure calculation failed.");
    }

#pragma omp parallel for
    for (auto& pi : particles) {
        pi.pressure = pressures[pi.id];
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
        if (pi.boundaryCondition == BoundaryCondition::Ignored)
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
            checkBoundaryViolation(pi);
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
        spdlog::debug(
            "Courant number is larger than CFL condition. Courant={}",
            std::to_string(maxCourantNumber)
        );
        debugLogCount++;
    }
    if (maxCourantNumber > 1.0) {
        spdlog::error("Courant number is larger than 1.0. Exiting Simulation.");
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

// Set the number density at the final position.
// This function should be called only in the following cases.
// 1. When the higher order source term for pressure calculation is on; this value will be
// used in the pressure calculation for the next time step.
// 2. When it's time to save; output the this value so that the correct number density can
// be displayed.
void MPS::setNumberDensity() {
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

void MPS::checkBoundaryViolation(Particle& pi) {
    bool isInDomain = true;
    if (pi.position.x() < settings.domain.x.min ||
        settings.domain.x.max < pi.position.x())
        isInDomain = false;
    if (pi.position.y() < settings.domain.y.min ||
        settings.domain.y.max < pi.position.y())
        isInDomain = false;
    if (pi.position.z() < settings.domain.z.min ||
        settings.domain.z.max < pi.position.z())
        isInDomain = false;
    if (!isInDomain) {
        spdlog::debug(
            "Particle (id={}) is out of domain. (x, y, z)=({}, {}, "
            "{}) ",
            std::to_string(pi.id),
            std::to_string(pi.position.x()),
            std::to_string(pi.position.y()),
            std::to_string(pi.position.z())
        );
        debugLogCount++;

        pi.type = ParticleType::Ghost;
    }
}
