#pragma once

#include <Eigen/Core>
#include <filesystem>
#include <fstream>

class Domain {
private:
public:
    struct {
        double min, max, length;
    } x, y, z;

    Domain() = default;
    Domain(
        const double& xMin,
        const double& xMax,
        const double& yMin,
        const double& yMax,
        const double& zMin,
        const double& zMax
    ) {
        this->x.min    = xMin;
        this->x.max    = xMax;
        this->x.length = xMax - xMin;
        this->y.min    = yMin;
        this->y.max    = yMax;
        this->y.length = yMax - yMin;
        this->z.min    = zMin;
        this->z.max    = zMax;
        this->z.length = zMax - zMin;
    }
};

class Settings {
private:
public:
    // computational conditions
    double dim;
    double particleDistance;
    double dt;
    double finishTime;
    double outputInterval;
    double cflCondition;

    // domain
    Domain domain;

    // physical properties
    double density;
    double kinematicViscosity;

    // gravity
    Eigen::Vector3d gravity;

    // effective radius
    struct {
        double pressure;
        double viscosity;
        double surfaceDetection;
        double max;
    } effectiveRadius;

    // surface detection
    double thresholdForSurfaceDetection;

    // pressure calculation
    struct {
        bool on;
        double compressibility;
    } quasiCompressibility;
    struct {
        bool on;
        double gamma;
    } higherOrderSourceTerm;
    double relaxationCoefficientForPressure;

    // collision
    double collisionDistance;
    double coefficientOfRestitution;

    // io
    std::filesystem::path inputCsvPath;

    void load(const std::filesystem::path& inputYamlPath);
};
