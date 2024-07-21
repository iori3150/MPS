#pragma once

#include <Eigen/Core>
#include <fstream>

struct Settings {
    // calculation conditions
    double dim                = 2;
    double particleDistance   = 0.05;
    double dt                 = 0.002;
    double outputInterval     = 0.05;
    double finishTime         = 2.0;
    double cflCondition       = 0.2;
    std::string inputProfPath = "input/result/input.prof";
    std::string inputDataPath = "input/result/input.data";

    // physical properties
    double kinematicViscosity = 1.0e-06;
    double density            = 1000.0;

    // effective radius
    struct {
        double numberDensity;
        double gradient;
        double laplacian;
    } re{2.1 * particleDistance, 2.1 * particleDistance, 3.1 * particleDistance};

    // gravity
    Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -9.8, 0.0);

    // collision
    double collisionDistance        = 0.5 * particleDistance;
    double coefficientOfRestitution = 0.2;

    // boundary condition
    double thresholdForSurfaceDetection = 0.97;

    // source term
    double relaxationCoefficientForPressure = 0.2;

    // matrix
    double compressibility = 0.45e-9;
};
