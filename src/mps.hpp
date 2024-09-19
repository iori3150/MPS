#pragma once

#include "bucket.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <Eigen/Sparse>
#include <map>

class RefValues {
private:
public:
    double initialNumberDensity = 0;
    double lambda               = 0;

    RefValues() = default;
    RefValues(
        const int& dim, const double& particleDistance, const double& effectiveRadius
    );
};

class MPS {
private:
    Bucket bucket;

    struct {
        RefValues pressure;
        RefValues viscosity;
        RefValues surfaceDetection;
    } refValues;

    Eigen::SparseMatrix<double, Eigen::RowMajor> coefficientMatrix;
    Eigen::VectorXd sourceTerm;

    double importInitialCondition(); // Import initial condition and return initial time

    void calcGravity();
    void calcViscosity();
    void moveParticles();
    void collision();

    void calcPressure();
    void setBoundaryCondition();
    bool isSurfaceParticle(Particle& pi);
    void setSourceTerm();
    void setMatrix();
    void ensureDirichletBoundaryConnection();
    void solvePoissonEquation();
    void removeNegativePressure();
    void setMinimumPressure();

    void calcPressureGradient();
    void moveParticlesWithPressureGradient();

    double getNumberDensity(const Particle& pi, const double& re);
    void setNumberDensity();
    void setNeighbors();
    void checkBoundaryViolation(Particle& pi);

    double getCourantNumber();

public:
    Settings settings;
    std::vector<Particle> particles;
    int debugLogCount = 0;

    MPS() = default;
    MPS(const std::filesystem::path& inputYamlPath);

    double loadInitialState(); // Initialize particles and return initial time
    double stepForward(const bool isTimeToExport);
};
