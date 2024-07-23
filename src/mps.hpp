#pragma once

#include "bucket.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <map>

class MPS {
private:
    struct {
        double numberDensity = 0.0;
        double gradient      = 0.0;
        double laplacian     = 0.0;
    } n0;
    double lambda = 0.0;

    Settings settings;
    Bucket bucket;

    Eigen::VectorXd sourceTerm;
    Eigen::MatrixXd coeffMatrix;
    std::vector<double> flagForCheckingBoundaryCondition;

    void calcGravity();
    void calcViscosity();
    void moveParticles();
    void collision();

    void calcPressure();
    void calcNumberDensity();
    void setBoundaryCondition();
    void setSourceTerm();
    void setMatrix();
    void exceptionalProcessingForBoundaryCondition();
    void checkBoundaryCondition();
    void increaseDiagonalTerm();
    void solvePoissonEquation();
    void removeNegativePressure();
    void setMinimumPressure();

    void calcPressureGradient();
    void moveParticlesWithPressureGradient();

    double weight(const double& dist, const double& re);

    void setNeighbors();

public:
    std::vector<Particle> particles;

    MPS() = default;
    MPS(const Settings& settings, std::vector<Particle>& particles);

    void stepForward();

    double getCourantNumber();
};
