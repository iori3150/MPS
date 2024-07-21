#pragma once

#include "mps.hpp"
#include "particle.hpp"

#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;

    void startSimulation();
    void endSimulation();

    // main()
    void read_data(std::vector<Particle>& particles);
    void set_parameter();
    void set_bucket();
    void main_loop();

    // main_loop()
    void write_data();

    // bucket
    void store_particle();
    void setNeighbors();

    // time calculation
    std::tuple<int, int, int> cal_h_m_s(int second);
};
