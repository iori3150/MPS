#include "simulation.hpp"

#include "bucket.hpp"
#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <omp.h>

using std::cerr;
using std::cout;
using std::endl;
namespace fs = std::filesystem;

#define rep(i, a, b) for (int i = a; i < b; i++)
#define ON 1
#define OFF 0

int np; // number of particles

double re_max, re2_max; // for bucket and neighor

// other parameters
double courantNumber;

// main()
clock_t sim_start_time;
int error_flag = OFF;

// main_loop()
int timestep;
double Time;
clock_t timestep_start_time;

// write_data()
int nfile; // number of files
FILE* log_file;

Settings settings;

void Simulation::run() {
    startSimulation();

    std::vector<Particle> particles;
    Domain domain;

    read_data(particles, domain);
    mps = MPS(settings, particles);

    set_parameter();
    bucket = Bucket(re_max, domain, particles.size());

    main_loop();

    endSimulation();
}

void Simulation::startSimulation() {
    cout << endl << "*** START SIMULATION ***" << endl;
    sim_start_time = clock();
}

void Simulation::endSimulation() {
    int hour, minute, second;
    second                         = (double) (clock() - sim_start_time) / CLOCKS_PER_SEC;
    std::tie(hour, minute, second) = cal_h_m_s(second);
    printf("\nTotal Simulation Time = %dh %02dm %02ds\n", hour, minute, second);

    if (error_flag == ON)
        cout << "Error has occured. Please check error.log" << endl;
    else
        cout << "There was no error." << endl;

    fclose(log_file);

    cout << endl << "*** END SIMULATION ***" << endl << endl;
}

void Simulation::read_data(std::vector<Particle>& particles, Domain& domain) {
    std::ifstream file;

    file.open(settings.inputProfPath);
    if (!file) {
        cout << "ERROR: There is no file named " << settings.inputProfPath << endl;
        exit(1);
    }

    file >> Time;
    file >> np;

    int id;
    rep(i, 0, np) {
        int type;
        double x, y, z, u, v, w;
        double pressure, n;

        file >> id >> type;
        file >> x >> y >> z;
        file >> u >> v >> w;
        file >> pressure >> n;

        particles.push_back(Particle(
            id,
            static_cast<ParticleType>(type),
            Eigen::Vector3d(x, y, z),
            Eigen::Vector3d(u, v, w),
            pressure,
            settings.density
        ));
    }

    file.close();

    file.open(settings.inputDataPath);
    if (!file) {
        cout << "ERROR: There is no file named " << settings.inputDataPath << endl;
        exit(1);
    }

    std::string dummy_string;
    double xMin, xMax, yMin, yMax, zMin, zMax;
    file >> dummy_string >> xMin;
    file >> dummy_string >> xMax;
    file >> dummy_string >> yMin;
    file >> dummy_string >> yMax;
    file >> dummy_string >> zMin;
    file >> dummy_string >> zMax;
    domain = Domain(xMin, xMax, yMin, yMax, zMin, zMax);

    file.close();
}

void Simulation::set_parameter() {

    // effective radius;
    re_max =
        std::max({settings.re.numberDensity, settings.re.gradient, settings.re.laplacian}
        );
    re2_max = re_max * re_max;

    // main_loop()
    Time = 0.0;

    // write_data()
    nfile = 0;
    char filename[256];
    sprintf(filename, "result/result.log");
    log_file = fopen(filename, "w");
}

void Simulation::main_loop() {
    timestep = 0;

    write_data();

    while (Time <= settings.finishTime) {
        timestep_start_time = clock();

        setNeighbors();
        mps.calGravity();
        mps.calViscosity();
        mps.moveParticle();

        setNeighbors();
        mps.collision();

        setNeighbors();
        mps.calcPressure();
        mps.calcPressureGradient();
        mps.moveParticleWithPressureGradient();

        courantNumber = mps.calcCourantNumber();

        timestep++;
        Time += settings.dt;
        write_data();
    }
}

void Simulation::write_data() {
    clock_t now = clock();
    int hour, minute, second;

    // elapsed
    char elapsed[256];
    second                         = (double) (now - sim_start_time) / CLOCKS_PER_SEC;
    std::tie(hour, minute, second) = cal_h_m_s(second);
    sprintf(elapsed, "elapsed=%dh %02dm %02ds", hour, minute, second);

    // ave [s]/[timestep]
    double ave = ((double) (now - sim_start_time) / CLOCKS_PER_SEC) / timestep;
    if (timestep == 0)
        ave = 0;

    // remain
    char remain[256];
    second = ((settings.finishTime - Time) / Time) * ave * timestep;
    std::tie(hour, minute, second) = cal_h_m_s(second);
    if (timestep == 0)
        sprintf(remain, "remain=---");
    else
        sprintf(remain, "remain=%dh %02dm %02ds", hour, minute, second);

    // last
    double last = (double) (now - timestep_start_time) / CLOCKS_PER_SEC;

    // terminal output
    printf(
        "%d: settings.dt=%.gs   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs/step   "
        "last=%.3lfs/step   out=%dfiles   Courant=%.2lf\n",
        timestep,
        settings.dt,
        Time,
        settings.finishTime,
        elapsed,
        remain,
        ave,
        last,
        nfile,
        courantNumber
    );

    // log file output
    fprintf(
        log_file,
        "%d: settings.dt=%gs   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs/step   "
        "last=%.3lfs/step   out=%dfiles   Courant=%.2lf\n",
        timestep,
        settings.dt,
        Time,
        settings.finishTime,
        elapsed,
        remain,
        ave,
        last,
        nfile,
        courantNumber
    );

    // error file output
    fprintf(stderr, "%4d: t=%.3lfs\n", timestep, Time);

    // prof file output
    if (Time >= settings.outputInterval * double(nfile)) {
        FILE* fp;
        char filename[256];

        sprintf(filename, "result/prof/output_%04d.prof", nfile);
        fp = fopen(filename, "w");
        fprintf(fp, "%lf\n", Time);
        fprintf(fp, "%d\n", np);
        rep(i, 0, np) {
            if (mps.particles[i].type == ParticleType::Ghost)
                continue;

            fprintf(
                fp,
                "%4d %2d % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 09.3lf "
                "% 08.3lf\n",
                i,
                mps.particles[i].type,
                mps.particles[i].position.x(),
                mps.particles[i].position.y(),
                mps.particles[i].position.z(),
                mps.particles[i].velocity[0],
                mps.particles[i].velocity[1],
                mps.particles[i].velocity[2],
                mps.particles[i].pressure,
                mps.particles[i].numberDensity
            );
        }
        fclose(fp);

        nfile++;
    }
}

void Simulation::setNeighbors() {
    bucket.storeParticles(mps.particles);

#pragma omp parallel for
    for (auto& pi : mps.particles) {
        if (pi.type == ParticleType::Ghost)
            continue;

        pi.neighbors.clear();

        int ix = int((pi.position.x() - bucket.domain.x.min) / bucket.length) + 1;
        int iy = int((pi.position.y() - bucket.domain.y.min) / bucket.length) + 1;
        int iz = int((pi.position.z() - bucket.domain.z.min) / bucket.length) + 1;

        for (int jx = ix - 1; jx <= ix + 1; jx++) {
            for (int jy = iy - 1; jy <= iy + 1; jy++) {
                for (int jz = iz - 1; jz <= iz + 1; jz++) {
                    int jBucket = jx + jy * bucket.numX + jz * bucket.numX * bucket.numY;
                    int j       = bucket.first[jBucket];

                    while (j != -1) {
                        Particle& pj = mps.particles[j];

                        double dist = (pj.position - pi.position).norm();
                        if (j != pi.id && dist < re_max) {
                            pi.neighbors.emplace_back(j, dist);
                        }

                        j = bucket.next[j];
                    }
                }
            }
        }
    }
}

std::tuple<int, int, int> Simulation::cal_h_m_s(int second) {
    int hour = second / 3600;
    second %= 3600;
    int minute = second / 60;
    second %= 60;

    return std::forward_as_tuple(hour, minute, second);
}
