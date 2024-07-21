#include "simulation.hpp"

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

// bucket
double x_min, x_max, y_min, y_max, z_min, z_max;
int num_bucket, num_bucket_x, num_bucket_y, num_bucket_xy, num_bucket_z;
double bucket_length;
Eigen::VectorXi bucket_next, bucket_first, bucket_last;

Settings settings;

void Simulation::run() {
    startSimulation();

    std::vector<Particle> particles;
    read_data(particles);
    set_parameter();
    set_bucket();

    main_loop(particles);

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

void read_data(std::vector<Particle>& particles) {
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
    file >> dummy_string >> x_min;
    file >> dummy_string >> x_max;
    file >> dummy_string >> y_min;
    file >> dummy_string >> y_max;
    file >> dummy_string >> z_min;
    file >> dummy_string >> z_max;

    file.close();
}

void set_parameter() {

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

void set_bucket() {
    bucket_length = re_max * (1.0 + settings.cflCondition);

    num_bucket_x  = (int) ((x_max - x_min) / bucket_length) + 3;
    num_bucket_y  = (int) ((y_max - y_min) / bucket_length) + 3;
    num_bucket_z  = (int) ((z_max - z_min) / bucket_length) + 3;
    num_bucket_xy = num_bucket_x * num_bucket_y;
    num_bucket    = num_bucket_x * num_bucket_y * num_bucket_z;

    bucket_first.resize(num_bucket);
    bucket_last.resize(num_bucket);
    bucket_next.resize(np);
}

void main_loop(std::vector<Particle>& particles) {
    timestep = 0;

    write_data(particles);

    MPS mps(settings, particles.size());

    while (Time <= settings.finishTime) {
        timestep_start_time = clock();

        store_particle(particles);
        setNeighbors(particles);
        mps.calGravity(particles);
        mps.calViscosity(particles);
        mps.moveParticle(particles);

        setNeighbors(particles);
        mps.collision(particles);

        setNeighbors(particles);
        mps.calcPressure(particles);
        mps.calcPressureGradient(particles);
        mps.moveParticleWithPressureGradient(particles);

        courantNumber = mps.calcCourantNumber(particles);

        timestep++;
        Time += settings.dt;
        write_data(particles);
    }
}

void write_data(std::vector<Particle>& particles) {
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
            if (particles[i].type == ParticleType::Ghost)
                continue;

            fprintf(
                fp,
                "%4d %2d % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 09.3lf "
                "% 08.3lf\n",
                i,
                particles[i].type,
                particles[i].position.x(),
                particles[i].position.y(),
                particles[i].position.z(),
                particles[i].velocity[0],
                particles[i].velocity[1],
                particles[i].velocity[2],
                particles[i].pressure,
                particles[i].numberDensity
            );
        }
        fclose(fp);

        nfile++;
    }
}

void store_particle(std::vector<Particle>& particles) {
#pragma omp parallel for
    rep(i, 0, num_bucket) {
        bucket_first[i] = -1;
        bucket_last[i]  = -1;
    }

#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].type == ParticleType::Ghost)
            continue;
        bucket_next[i] = -1;
    }

#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].type == ParticleType::Ghost)
            continue;

        int ix      = (int) ((particles[i].position.x() - x_min) / bucket_length) + 1;
        int iy      = (int) ((particles[i].position.y() - y_min) / bucket_length) + 1;
        int iz      = (int) ((particles[i].position.z() - z_min) / bucket_length) + 1;
        int ibucket = iz * num_bucket_xy + iy * num_bucket_x + ix;

#pragma omp critical
        {
            if (bucket_last[ibucket] == -1)
                bucket_first[ibucket] = i;
            else
                bucket_next[bucket_last[ibucket]] = i;
            bucket_last[ibucket] = i;
        }
    }
}

void setNeighbors(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost)
            continue;

        pi.neighbors.clear();

        int ix = int((pi.position.x() - x_min) / bucket_length) + 1;
        int iy = int((pi.position.y() - y_min) / bucket_length) + 1;
        int iz = int((pi.position.z() - z_min) / bucket_length) + 1;

        for (int jx = ix - 1; jx <= ix + 1; jx++) {
            for (int jy = iy - 1; jy <= iy + 1; jy++) {
                for (int jz = iz - 1; jz <= iz + 1; jz++) {
                    int jbucket = jz * num_bucket_xy + jy * num_bucket_x + jx;
                    int j       = bucket_first[jbucket];

                    while (j != -1) {
                        Particle& pj = particles[j];

                        double dist = (pj.position - pi.position).norm();
                        if (j != pi.id && dist < re_max) {
                            pi.neighbors.emplace_back(j, dist);
                        }

                        j = bucket_next[j];
                    }
                }
            }
        }
    }
}

std::tuple<int, int, int> cal_h_m_s(int second) {
    int hour = second / 3600;
    second %= 3600;
    int minute = second / 60;
    second %= 60;

    return std::forward_as_tuple(hour, minute, second);
}
