#include "bucket.hpp"

#include "particle.hpp"

#include <iostream>

using std::cerr;
using std::endl;

Bucket::Bucket(const double& reMax, const Domain& domain, const int& particleSize) {
    this->length = reMax;
    this->domain = domain;

    this->numX = (int) (domain.x.length / length) + 3;
    this->numY = (int) (domain.y.length / length) + 3;
    this->numZ = (int) (domain.z.length / length) + 3;
    this->num  = numX * numY * numZ;

    this->first.resize(num);
    this->last.resize(num);
    this->next.resize(particleSize);
}

void Bucket::storeParticles(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (int i = 0; i < num; i++) {
        first[i] = -1;
        last[i]  = -1;
    }

#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        next[i] = -1;
    }

    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost)
            continue;

        int ix      = (int) ((pi.position.x() - domain.x.min) / length) + 1;
        int iy      = (int) ((pi.position.y() - domain.y.min) / length) + 1;
        int iz      = (int) ((pi.position.z() - domain.z.min) / length) + 1;
        int iBucket = ix + iy * numX + iz * numX * numY;

        if (last[iBucket] == -1)
            first[iBucket] = pi.id;
        else
            next[last[iBucket]] = pi.id;
        last[iBucket] = pi.id;
    }
}
