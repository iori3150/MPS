#pragma once

#include "mps.hpp"
#include "particle.hpp"
#include "saver.hpp"
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
    std::chrono::system_clock::time_point simulationStartTime;

    int timestep = 0;
    double time  = 0.0;

    int resultFileNum = 0;
    std::ofstream logFile;
    Saver saver;

    void startSimulation();
    void endSimulation();

    void read_data(std::vector<Particle>& particles);
    void timeStepReport(
        const std::chrono::system_clock::time_point& timeStepStartTime,
        const std::chrono::system_clock::time_point& timeStepEndTime,
        const double& courantNumber
    );
};
