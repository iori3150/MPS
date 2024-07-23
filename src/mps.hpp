#pragma once

#include "bucket.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <map>

class MPS {
private:
    struct {
        double pressure         = 0.0;
        double viscosity        = 0.0;
        double surfaceDetection = 0.0;
    } initialNumberDensity;
    struct {
        double pressure         = 0.0;
        double viscosity        = 0.0;
        double surfaceDetection = 0.0;
    } lambda;

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

    double getNumberDensity(const Particle& pi, const double& re);
    double weight(const double& dist, const double& re);

    void setNeighbors();

    void setNumberDensityForDisplay();

public:
    std::vector<Particle> particles;

    MPS() = default;
    MPS(const Settings& settings, std::vector<Particle>& particles);

    void stepForward(const bool isTimeToSave);

    double getCourantNumber();
};
