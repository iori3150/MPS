#pragma once

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

    Eigen::VectorXd sourceTerm;
    Eigen::MatrixXd coeffMatrix;
    std::vector<double> flagForCheckingBoundaryCondition;

    double weight(const double& dist, const double& re);

public:
    std::vector<Particle> particles;

    MPS(){};
    MPS(const Settings& settings, std::vector<Particle>& particles);

    void calGravity();
    void calViscosity();
    void moveParticle();
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
    void moveParticleWithPressureGradient();

    double calcCourantNumber();
};
