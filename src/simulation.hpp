#pragma once

#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <chrono>
#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;
    Settings settings;
    std::chrono::system_clock::time_point startTime;

    int timestep = 0;
    double time  = 0.0;

    int resultFileNum;
    FILE* logFile;

    void startSimulation();
    void endSimulation();

    // main()
    void read_data(std::vector<Particle>& particles);
    void set_parameter();

    // main_loop()
    void write_data(
        const double& courantNumber,
        const std::chrono::system_clock::time_point& timestepStartTime
    );

    std::tuple<int, int, int> getTimeDuration(
        const std::chrono::system_clock::time_point& start,
        const std::chrono::system_clock::time_point& end
    );
};
