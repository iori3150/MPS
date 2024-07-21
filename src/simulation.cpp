#include "simulation.hpp"

#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
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

// check_boundary_condition()
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0
#define DIRICHLET_BOUNDARY_IS_CONNECTED 1
#define DIRICHLET_BOUNDARY_IS_CHECKED 2

// main()
void read_data();
void set_parameter();
void cal_n0_and_lambda();
void set_bucket();
void main_loop();

// main_loop()
void write_data();
void cal_P(MPS& mps);
void cal_P_grad();
void move_particle_using_P_grad();
void cal_courant();

// cal_P()
void set_matrix();
void exceptional_processing_for_boundary_condition();
void check_boundary_condition();
void increase_diagonal_term();
void solve_Poisson_eq(MPS& mps);
void remove_negative_P();
void set_P_min();

// bucket
void store_particle();
void setNeighbors(std::vector<Particle>& particles);

// common
double weight(double dis, double re);
double cal_dis2(int i, int j);

// time calculation
std::tuple<int, int, int> cal_h_m_s(int second);

// particles
std::vector<Particle> particles;

int np; // number of particles
Eigen::VectorXd P_min;

// effective radius
double re_for_n, re2_for_n;
double re_for_grad, re2_for_grad;
double re_for_lap, re2_for_lap;
double re_max, re2_max; // for bucket and neighor

// other parameters
double n0_for_n;
double n0_for_grad;
double n0_for_lap;
double lambda; // used for laplacian calculation
double courant;

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

// cal_P()
Eigen::VectorXi flag_for_checking_boundary_condition;
Eigen::MatrixXd coef_matrix;

// bucket
double x_min, x_max, y_min, y_max, z_min, z_max;
int num_bucket, num_bucket_x, num_bucket_y, num_bucket_xy, num_bucket_z;
double bucket_length;
Eigen::VectorXi bucket_next, bucket_first, bucket_last;

Settings settings;

void Simulation::run() {
    startSimulation();

    read_data();
    set_parameter();
    set_bucket();

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

void read_data() {
    std::ifstream file;

    file.open(settings.inputProfPath);
    if (!file) {
        cout << "ERROR: There is no file named " << settings.inputProfPath << endl;
        exit(1);
    }

    file >> Time;
    file >> np;

    // set std::vector size
    P_min.resize(np);

    flag_for_checking_boundary_condition.resize(np);
    coef_matrix.resize(np, np);

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
    re_for_n     = settings.re.numberDensity;
    re2_for_n    = re_for_n * re_for_n;
    re_for_grad  = settings.re.gradient;
    re2_for_grad = re_for_grad * re_for_grad;
    re_for_lap   = settings.re.laplacian;
    re2_for_lap  = re_for_lap * re_for_lap;
    re_max       = std::max({re_for_n, re_for_grad, re_for_lap});
    re2_max      = re_max * re_max;

    // main_loop()
    Time = 0.0;

    // write_data()
    nfile = 0;
    char filename[256];
    sprintf(filename, "result/result.log");
    log_file = fopen(filename, "w");

    cal_n0_and_lambda();
}

void cal_n0_and_lambda() {
    n0_for_n    = 0.0;
    n0_for_grad = 0.0;
    n0_for_lap  = 0.0;
    lambda      = 0.0;

    int iZ_start, iZ_end;
    if (settings.dim == 2) {
        iZ_start = 0;
        iZ_end   = 1;
    } else {
        iZ_start = -4;
        iZ_end   = 5;
    }

    double xi, yi, zi;
    double dis, dis2;
    for (int iX = -4; iX < 5; iX++) {
        for (int iY = -4; iY < 5; iY++) {
            for (int iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))
                    continue;

                xi   = settings.particleDistance * (double) (iX);
                yi   = settings.particleDistance * (double) (iY);
                zi   = settings.particleDistance * (double) (iZ);
                dis2 = xi * xi + yi * yi + zi * zi;

                dis = sqrt(dis2);

                n0_for_n += weight(dis, re_for_n);
                n0_for_grad += weight(dis, re_for_grad);
                n0_for_lap += weight(dis, re_for_lap);

                lambda += dis2 * weight(dis, re_for_lap);
            }
        }
    }
    lambda /= n0_for_lap;
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

void main_loop() {
    timestep = 0;

    write_data();

    MPS mps(settings, particles.size());

    while (Time <= settings.finishTime) {
        timestep_start_time = clock();

        // explicit
        store_particle();
        setNeighbors(particles);
        mps.calGravity(particles);
        mps.calViscosity(particles);
        mps.moveParticle(particles);

        setNeighbors(particles);
        mps.collision(particles);

        // inplicit
        setNeighbors(particles);
        cal_P(mps);
        cal_P_grad();
        move_particle_using_P_grad();

        cal_courant();

        timestep++;
        Time += settings.dt;
        write_data();
    }
}

void write_data() {
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
        courant
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
        courant
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

void cal_P(MPS& mps) {
    mps.calcNumberDensity(particles);
    mps.setBoundaryCondition(particles);
    mps.setSourceTerm(particles);
    set_matrix();
    solve_Poisson_eq(mps);
    remove_negative_P();
    set_P_min();
}

void set_matrix() {
    coef_matrix.setZero();

    double n0 = n0_for_lap;
    double A  = 2.0 * settings.dim / (n0 * lambda);
#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].boundaryCondition != BoundaryCondition::Inner)
            continue;

        for (auto& neighbor : particles[i].neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re_for_lap) {
                double coef_ij = A * weight(dist, re_for_lap) / particles[i].density;

                coef_matrix(i, pj.id) = (-1.0 * coef_ij);
                coef_matrix(i, i) += coef_ij;
            }
        }

        coef_matrix(i, i) += (settings.compressibility) / (settings.dt * settings.dt);
    }

    exceptional_processing_for_boundary_condition();
}

void exceptional_processing_for_boundary_condition() {
    // If tere is no Dirichlet boundary condition on the fluid,
    // increase the diagonal terms of the matrix for an exception.
    // This allows us to solve the matrix without Dirichlet boundary conditions.
    check_boundary_condition();
    increase_diagonal_term();
}

void check_boundary_condition() {
#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].boundaryCondition == BoundaryCondition::GhostOrDummy) {
            flag_for_checking_boundary_condition[i] = -1;

        } else if (particles[i].boundaryCondition == BoundaryCondition::Surface) {
            flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_CONNECTED;

        } else {
            flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
        }
    }

    int count;
    while (true) {
        count = 0;

        rep(i, 0, np) {
            if (flag_for_checking_boundary_condition[i] !=
                DIRICHLET_BOUNDARY_IS_CONNECTED)
                continue;

            for (auto& neighbor : particles[i].neighbors) {
                const Particle& pj = particles[neighbor.id];
                const double& dist = neighbor.distance;

                if (flag_for_checking_boundary_condition[pj.id] !=
                    DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                    continue;

                if (dist < re_for_lap)
                    flag_for_checking_boundary_condition[pj.id] =
                        DIRICHLET_BOUNDARY_IS_CONNECTED;
            }

            flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
            count++;
        }

        if (count == 0)
            break;
    }

#pragma omp parallel for
    rep(i, 0, np) {
        if (flag_for_checking_boundary_condition[i] ==
            DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
#pragma omp critical
            {
                cerr << "WARNING: There is no dirichlet boundary condition for particle "
                     << i << endl;
            }
        }
    }
}

void increase_diagonal_term() {
#pragma omp parallel for
    rep(i, 0, np) {
        if (flag_for_checking_boundary_condition[i] ==
            DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
            coef_matrix(i, i) *= 2.0;
        }
    }
}

void solve_Poisson_eq(MPS& mps) {
    Eigen::SparseMatrix<double> A = coef_matrix.sparseView();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd pressures = solver.solve(mps.sourceTerm);
    for (auto& pi : particles) {
        pi.pressure = pressures[pi.id];
    }
}

void remove_negative_P() {
#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].pressure < 0.0)
            particles[i].pressure = 0.0;
    }
}

void set_P_min() {
#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].boundaryCondition == BoundaryCondition::GhostOrDummy)
            continue;

        P_min[i] = particles[i].pressure;

        for (auto& neighbor : particles[i].neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re_for_grad) {
                if (P_min[i] > pj.pressure)
                    P_min[i] = pj.pressure;
            }
        }
    }
}

void cal_P_grad() {
    double A = settings.dim / n0_for_grad;

#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].type != ParticleType::Fluid)
            continue;

        Eigen::Vector3d grad = Eigen::Vector3d::Zero();

        for (auto& neighbor : particles[i].neighbors) {
            const Particle& pj = particles[neighbor.id];
            const double& dist = neighbor.distance;

            if (pj.type == ParticleType::DummyWall)
                continue;

            if (dist < re_for_grad) {
                grad += (pj.position - particles[i].position) * (pj.pressure - P_min[i]) *
                        weight(dist, re_for_grad) / (dist * dist);
            }
        }

        grad *= A;
        particles[i].acceleration = -1.0 * grad / particles[i].density;
    }
}

void move_particle_using_P_grad() {
#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].type == ParticleType::Fluid) {
            particles[i].velocity += particles[i].acceleration * settings.dt;
            particles[i].position +=
                particles[i].acceleration * settings.dt * settings.dt;
        }

        particles[i].acceleration.Zero();
    }
}

void cal_courant() {
    courant = 0.0;

#pragma omp parallel for
    rep(i, 0, np) {
        if (particles[i].type != ParticleType::Fluid)
            continue;

        double courant_i =
            (particles[i].velocity.norm() * settings.dt) / settings.particleDistance;
        if (courant_i > courant)
#pragma omp critical
        {
            courant = courant_i;
        }
    }

    if (courant > settings.cflCondition) {
        cerr << "ERROR: Courant number is larger than CFL condition. Courant = "
             << courant << endl;
        error_flag = ON;
    }
}

void store_particle() {
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

double weight(double dis, double re) {
    double w = 0.0;

    if (dis < re)
        w = (re / dis) - 1.0;

    return w;
}

double cal_dis2(int i, int j) {
    Eigen::Vector3d x_ij = particles[j].position - particles[i].position;
    return x_ij.squaredNorm();
}

std::tuple<int, int, int> cal_h_m_s(int second) {
    int hour = second / 3600;
    second %= 3600;
    int minute = second / 60;
    second %= 60;

    return std::forward_as_tuple(hour, minute, second);
}
