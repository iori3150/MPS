#pragma once

#include "exporter.hpp"
#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <chrono>
#include <fstream>
#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;
    Settings settings;

    int timeStep = 0;
    double time  = 0.0;
    std::chrono::system_clock::time_point simulationStartTime, simulationEndTime;

    Exporter exporter;
    int outFileNum = 0;

    std::ofstream logFile;

    void startSimulation();
    void endSimulation();

    void read_data(std::vector<Particle>& particles);
    void timeStepReport(
        const std::chrono::system_clock::time_point& timeStepStartTime,
        const std::chrono::system_clock::time_point& timeStepEndTime,
        const double& courantNumber
    );
    bool isTimeToExport();
    void exportParticles(const std::vector<Particle>& particles);
};
