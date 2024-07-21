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

    std::vector<double> flagForCheckingBoundaryCondition;

    double weight(const double& dist, const double& re);

public:
    Eigen::VectorXd sourceTerm;
    Eigen::MatrixXd coeffMatrix;

    MPS(const Settings& settings, const int& numberOfParticles);

    void calGravity(std::vector<Particle>& particles);
    void calViscosity(std::vector<Particle>& particles);
    void moveParticle(std::vector<Particle>& particles);
    void collision(std::vector<Particle>& particles);

    void calcNumberDensity(std::vector<Particle>& particles);
    void setBoundaryCondition(std::vector<Particle>& particles);
    void setSourceTerm(std::vector<Particle>& particles);
    void setMatrix(std::vector<Particle>& particles);
    void exceptionalProcessingForBoundaryCondition(std::vector<Particle>& particles);
    void checkBoundaryCondition(std::vector<Particle>& particles);
    void increaseDiagonalTerm(std::vector<Particle>& particles);
};
