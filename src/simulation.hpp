#pragma once

#include "bucket.hpp"
#include "mps.hpp"
#include "particle.hpp"

#include <vector>

class Simulation {
public:
    void run();

private:
    MPS mps;
    Bucket bucket;

    void startSimulation();
    void endSimulation();

    // main()
    void read_data(std::vector<Particle>& particles);
    void set_parameter();
    void main_loop();

    // main_loop()
    void write_data();

    // bucket
    void setNeighbors();

    // time calculation
    std::tuple<int, int, int> cal_h_m_s(int second);
};
