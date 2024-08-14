#include "../../src/exporter.hpp"
#include "../../src/particle.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

double l0      = 0.005;
double density = 1000.0;
std::vector<double> domainRangeX{0.0, 0.2};
std::vector<double> domainRangeY{0.0, 0.35};
std::vector<double> fluidRangeX{0.0, 0.2};
std::vector<double> fluidRangeY{0.0, 0.3};
int wallLayer      = 3;
int dummyWallLayer = 1;

void checkFluidRange(std::vector<double>& xRange, std::vector<double>& yRange);
bool isInside(
    Eigen::Vector3d& position, std::vector<double>& xRange, std::vector<double>& yRange
);

int main() {
    std::vector<Particle> particles;

    checkFluidRange(fluidRangeX, fluidRangeY);

    int allWallLayer = dummyWallLayer + wallLayer;

    int ixBegin = round(domainRangeX[0] / l0) - allWallLayer;
    int ixEnd   = round(domainRangeX[1] / l0) + allWallLayer;
    int iyBegin = round(domainRangeY[0] / l0) - allWallLayer;
    int iyEnd   = round(domainRangeY[1] / l0) + allWallLayer;

    for (int ix = ixBegin; ix <= ixEnd; ix++) {
        for (int iy = iyBegin; iy <= iyEnd; iy++) {
            Eigen::Vector3d pos((ix + 0.5) * l0, (iy + 0.5) * l0, 0.0);
            ParticleType type = ParticleType::Ghost;
            std::vector<double> xRange(2), yRange(2);
            double xMax = domainRangeX[1];
            double yMax = domainRangeY[1];

            // dummy wall region
            xRange = {-allWallLayer * l0, xMax + allWallLayer * l0};
            yRange = {-allWallLayer * l0, yMax};
            if (isInside(pos, xRange, yRange))
                type = ParticleType::DummyWall;

            // wall region
            xRange = {-wallLayer * l0, xMax + wallLayer * l0};
            yRange = {-wallLayer * l0, yMax};
            if (isInside(pos, xRange, yRange))
                type = ParticleType::Wall;

            // empty region
            xRange = {0.0, xMax};
            yRange = {0.0, yMax};
            if (isInside(pos, xRange, yRange))
                type = ParticleType::Ghost;

            // fluid region
            if (isInside(pos, fluidRangeX, fluidRangeY))
                type = ParticleType::Fluid;

            if (type != ParticleType::Ghost) {
                Eigen::Vector3d vel = Eigen::Vector3d::Zero();
                particles.push_back(Particle(particles.size(), type, pos, vel, density));
            }
        }
    }

    Exporter exporter;
    exporter.toCsv(fs::path("input.csv"), particles, 0.0);
    exporter.toVtu(fs::path("input.vtu"), particles, 0.0);

    std::cout << "Input file was created successfully." << std::endl;
}

void checkFluidRange(std::vector<double>& xRange, std::vector<double>& yRange) {
    double eps = 0.01 * l0;

    std::sort(xRange.begin(), xRange.end());
    std::sort(yRange.begin(), yRange.end());

    double nx = (xRange[1] - xRange[0]) / l0;
    if (abs(nx - std::round(nx)) > eps) {
        std::cout << "ERROR: xRange of the fluid is not divisible by particle length."
                  << std::endl;
        std::exit(-1);
    }

    double ny = (yRange[1] - yRange[0]) / l0;
    if (abs(ny - std::round(ny)) > eps) {
        std::cout << "ERROR: yRange of the fluid is not divisible by particle length."
                  << std::endl;
        std::exit(-1);
    }
}

bool isInside(
    Eigen::Vector3d& pos, std::vector<double>& xRange, std::vector<double>& yRange
) {
    double eps = 0.01 * l0;

    std::sort(xRange.begin(), xRange.end());
    std::sort(yRange.begin(), yRange.end());

    return (xRange[0] - eps < pos.x() && pos.x() < xRange[1] + eps) &&
           (yRange[0] - eps < pos.y() && pos.y() < yRange[1] + eps);
}
