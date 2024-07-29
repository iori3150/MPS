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
    Settings settings;
    Bucket bucket;

    struct {
        RefValues pressure;
        RefValues viscosity;
        RefValues surfaceDetection;
    } refValues;

    Eigen::SparseMatrix<double, Eigen::RowMajor> coefficientMatrix;
    Eigen::VectorXd sourceTerm;

    void importInitialCondition(double& time);

    void calcGravity();
    void calcViscosity();
    void moveParticles();
    void collision();

    void calcPressure();
    void setBoundaryCondition();
    void setSourceTerm();
    void setMatrix();
    void ensureDirichletBoundaryConnection();
    void solvePoissonEquation();
    void removeNegativePressure();
    void setMinimumPressure();

    void calcPressureGradient();
    void moveParticlesWithPressureGradient();

    double getNumberDensity(const Particle& pi, const double& re);

    void setNeighbors();

    void setNumberDensityForDisplay();

public:
    std::vector<Particle> particles;

    MPS() = default;
    MPS(const Settings& settings, double& time);

    void stepForward(const bool isTimeToExport);

    double getCourantNumber();
};
