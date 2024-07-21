#pragma once

#include "particle.hpp"

#include <vector>

class Simulation {
public:
    void run();

private:
    void startSimulation();
    void endSimulation();

    // main()
    void read_data(std::vector<Particle>& particles);
    void set_parameter();
    void set_bucket();
    void main_loop(std::vector<Particle>& particles);

    // main_loop()
    void write_data(std::vector<Particle>& particles);

    // bucket
    void store_particle(std::vector<Particle>& particles);
    void setNeighbors(std::vector<Particle>& particles);

    // time calculation
    std::tuple<int, int, int> cal_h_m_s(int second);
};
