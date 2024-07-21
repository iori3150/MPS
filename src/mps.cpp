#include "mps.hpp"

#define rep(i, a, b) for (int i = a; i < b; i++)

MPS::MPS(Settings& settings, std::vector<Particle>& particles) {
    this->settings  = settings;
    this->particles = particles;

    int iZ_start = -4;
    int iZ_end   = 5;
    if (settings.dim == 2) {
        iZ_start = 0;
        iZ_end   = 1;
    }
    for (int iX = -4; iX < 5; iX++) {
        for (int iY = -4; iY < 5; iY++) {
            for (int iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))
                    continue;

                Eigen::Vector3d r(
                    settings.particleDistance * (double) iX,
                    settings.particleDistance * (double) iY,
                    settings.particleDistance * (double) iZ
                );
                double dist  = r.norm();
                double dist2 = dist * dist;

                n0["gradient"] += weight(dist, settings.re_forGradient);
                n0["laplacian"] += weight(dist, settings.re_forLaplacian);
                lambda += dist2 * weight(dist, settings.re_forLaplacian);
            }
        }
    }
    lambda /= n0["laplacian"];
}

void MPS::calGravity(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Fluid) {
            pi.acceleration = settings.gravity;

        } else {
            pi.acceleration.setZero();
        }
    }
}

// void MPS::calViscosity(std::vector<Particle>& particles) {
//     double A =
//         (settings.kinematicViscosity) * (2.0 * settings.dim) / (n0["lap"] * lambda);
//
// #pragma omp parallel for
//     for (auto& pi : particles) {
//         if (pi.type != ParticleType::Fluid)
//             continue;
//
//         Eigen::Vector3d viscosity_term = Eigen::Vector3d::Zero();
//
//         rep(j_neighbor, 0, num_neighbor[i]) {
//             int j       = neighbor_id[i][j_neighbor];
//             double dis2 = neighbor_dis2[i][j_neighbor];
//
//             if (dis2 < re2_for_lap) {
//                 double dis = sqrt(dis2);
//                 viscosity_term += (particles[j].velocity - particles[i].velocity) *
//                                   weight(dis, re_for_lap);
//             }
//         }
//
//         viscosity_term *= A;
//         pi.acceleration += viscosity_term;
//     }
// }

double MPS::weight(const double& dist, const double& re) {
    double w = 0.0;

    if (dist < re)
        w = (re / dist) - 1.0;

    return w;
}
