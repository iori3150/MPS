#pragma once

#include <Eigen/Core>
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

class Settings{
private:
  // computational conditions
    double dim;
    double particleDistance;
    double dt; 
    double finishTime;
    double outputPeriod;
    double cflCondition;
    std::string inputProfPath                  = "input/input.csv";
    std::string inputDataPath                  = "input/result/input.data";

// domain
  Domain domain;

    // physical properties
    double density      ;
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
  double compressibility;
  double relaxationCoefficientForPressure;

  // collision
  double collisionDistance;
  double coefficientOfRestitution;

public:
    void load();
};
