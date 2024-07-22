#pragma once

#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;
    Settings settings;

    void startSimulation();
    void endSimulation();

    // main()
    void read_data(std::vector<Particle>& particles);
    void set_parameter();

    // main_loop()
    void write_data(const double& courantNumber);

    // time calculation
    std::tuple<int, int, int> cal_h_m_s(int second);
};
