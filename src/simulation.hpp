#pragma once

#include "mps.hpp"
#include "particle.hpp"

#include <chrono>
#include <fstream>
#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;

    int timeStep = 0;
    double time  = 0.0;
    std::chrono::system_clock::time_point simulationStartTime, simulationEndTime;

    int outFileNum = 0;
    std::ofstream logFile;
    std::filesystem::path resultDirectory;

    void startSimulation();
    void endSimulation();

    void createResultDirectory();

    void timeStepReport(
        const std::chrono::system_clock::time_point& timeStepStartTime,
        const std::chrono::system_clock::time_point& timeStepEndTime,
        const double& courantNumber
    );
    bool isTimeToExport();
    void exportParticles(const std::vector<Particle>& particles);
};
