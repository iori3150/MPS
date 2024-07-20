#include "simulation.hpp"

#include <omp.h>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>
#include <iostream>

#include "particle.hpp"
#include "settings.hpp"

using std::cerr;
using std::cout;
using std::endl;
namespace fs = std::filesystem;

#define rep(i, a, b) for (int i = a; i < b; i++)
#define ON 1
#define OFF 0

// set_boundary_condition
#define GHOST_OR_DUMMY -1
#define SURFACE_PARTICLE 1
#define INNER_PARTICLE 0

// check_boundary_condition()
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0
#define DIRICHLET_BOUNDARY_IS_CONNECTED 1
#define DIRICHLET_BOUNDARY_IS_CHECKED 2

// neighbor
#define NEIGHBOR_ARRAY_SIZE 50

// main()
void read_data();
void set_parameter();
void cal_n0_and_lambda();
void set_bucket();
void main_loop();

// main_loop()
void write_data();
void cal_gravity();
void cal_viscosity();
void move_particle();
void collision();
void cal_P();
void cal_P_grad();
void move_particle_using_P_grad();
void cal_courant();

// cal_P()
void cal_n();
void set_boundary_condition();
void set_source_term();
void set_matrix();
void exceptional_processing_for_boundary_condition();
void check_boundary_condition();
void increase_diagonal_term();
void solve_Poisson_eq();
void remove_negative_P();
void set_P_min();

// bucket
void store_particle();
void set_neighbor();

// common
double weight(double dis, double re);
double cal_dis2(int i, int j);

// time calculation
std::tuple<int, int, int> cal_h_m_s(int second);

// particles
std::vector<Particle> particles;

int                          np;  // number of particles
std::vector<Eigen::Vector3d> a;
Eigen::VectorXd              P, P_min, n;

// effective radius
double re_for_n, re2_for_n;
double re_for_grad, re2_for_grad;
double re_for_lap, re2_for_lap;
double re_max, re2_max;  // for bucket and neighor

// other parameters
double rho;
double n0_for_n;
double n0_for_grad;
double n0_for_lap;
double lambda;  // used for laplacian calculation
double courant;

// main()
clock_t sim_start_time;
int     error_flag = OFF;

// main_loop()
int     timestep;
double  Time;
clock_t timestep_start_time;

// write_data()
int   nfile;  // number of files
FILE *log_file;

// collision()
double collision_dis, collision_dis2;

// cal_P()
Eigen::VectorXi boundary_condition, flag_for_checking_boundary_condition;
Eigen::VectorXd source_term;
Eigen::MatrixXd coef_matrix;

// bucket
double          x_min, x_max, y_min, y_max, z_min, z_max;
int             num_bucket, num_bucket_x, num_bucket_y, num_bucket_xy, num_bucket_z;
double          bucket_length;
Eigen::VectorXi bucket_next, bucket_first, bucket_last;

// neighbor
Eigen::VectorXi              num_neighbor;
std::vector<Eigen::VectorXi> neighbor_id;
std::vector<Eigen::VectorXd> neighbor_dis2;

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
  cout << endl
       << "*** START SIMULATION ***" << endl;
  sim_start_time = clock();
}

void Simulation::endSimulation() {
  int hour, minute, second;
  second                         = (double)(clock() - sim_start_time) / CLOCKS_PER_SEC;
  std::tie(hour, minute, second) = cal_h_m_s(second);
  printf("\nTotal Simulation Time = %dh %02dm %02ds\n", hour, minute, second);

  if (error_flag == ON)
    cout << "Error has occured. Please check error.log" << endl;
  else
    cout << "There was no error." << endl;

  fclose(log_file);

  cout << endl
       << "*** END SIMULATION ***" << endl
       << endl;
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
  a.resize(np);
  P.resize(np);
  P_min.resize(np);
  n.resize(np);

  num_neighbor.resize(np);
  neighbor_id.resize(np);
  rep(i, 0, np) neighbor_id[i].resize(NEIGHBOR_ARRAY_SIZE);
  neighbor_dis2.resize(np);
  rep(i, 0, np) neighbor_dis2[i].resize(NEIGHBOR_ARRAY_SIZE);

  boundary_condition.resize(np);
  source_term.resize(np);
  flag_for_checking_boundary_condition.resize(np);
  coef_matrix.resize(np, np);

  int id;
  rep(i, 0, np) {
    int    type;
    double x, y, z, u, v, w;
    file >> id >> type;
    file >> x >> y >> z;
    file >> u >> v >> w;
    file >> P[i] >> n[i];
    particles.push_back(Particle(id, static_cast<ParticleType>(type), Eigen::Vector3d(x, y, z), Eigen::Vector3d(u, v, w)));
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
  rho = settings.density;

  // effective radius;
  re_for_n     = settings.re_forNumberDensity;
  re2_for_n    = re_for_n * re_for_n;
  re_for_grad  = settings.re_forGradient;
  re2_for_grad = re_for_grad * re_for_grad;
  re_for_lap   = settings.re_forLaplacian;
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

  // collision()
  collision_dis  = settings.collisionDistance;
  collision_dis2 = collision_dis * collision_dis;

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
        if (((iX == 0) && (iY == 0)) && (iZ == 0)) continue;

        xi   = settings.particleDistance * (double)(iX);
        yi   = settings.particleDistance * (double)(iY);
        zi   = settings.particleDistance * (double)(iZ);
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

  num_bucket_x  = (int)((x_max - x_min) / bucket_length) + 3;
  num_bucket_y  = (int)((y_max - y_min) / bucket_length) + 3;
  num_bucket_z  = (int)((z_max - z_min) / bucket_length) + 3;
  num_bucket_xy = num_bucket_x * num_bucket_y;
  num_bucket    = num_bucket_x * num_bucket_y * num_bucket_z;

  bucket_first.resize(num_bucket);
  bucket_last.resize(num_bucket);
  bucket_next.resize(np);
}

void main_loop() {
  timestep = 0;

  write_data();

  while (Time <= settings.finishTime) {
    timestep_start_time = clock();

    // explicit
    store_particle();
    set_neighbor();
    cal_gravity();
    cal_viscosity();
    move_particle();

    set_neighbor();
    collision();

    // inplicit
    set_neighbor();
    cal_P();
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
  int     hour, minute, second;

  // elapsed
  char elapsed[256];
  second                         = (double)(now - sim_start_time) / CLOCKS_PER_SEC;
  std::tie(hour, minute, second) = cal_h_m_s(second);
  sprintf(elapsed, "elapsed=%dh %02dm %02ds", hour, minute, second);

  // ave [s]/[timestep]
  double ave = ((double)(now - sim_start_time) / CLOCKS_PER_SEC) / timestep;
  if (timestep == 0) ave = 0;

  // remain
  char remain[256];
  second                         = ((settings.finishTime - Time) / Time) * ave * timestep;
  std::tie(hour, minute, second) = cal_h_m_s(second);
  if (timestep == 0)
    sprintf(remain, "remain=---");
  else
    sprintf(remain, "remain=%dh %02dm %02ds", hour, minute, second);

  // last
  double last = (double)(now - timestep_start_time) / CLOCKS_PER_SEC;

  // terminal output
  printf("%d: settings.dt=%.gs   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs/step   last=%.3lfs/step   out=%dfiles   Courant=%.2lf\n",
         timestep, settings.dt, Time, settings.finishTime, elapsed, remain, ave, last, nfile, courant);

  // log file output
  fprintf(log_file, "%d: settings.dt=%gs   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs/step   last=%.3lfs/step   out=%dfiles   Courant=%.2lf\n",
          timestep, settings.dt, Time, settings.finishTime, elapsed, remain, ave, last, nfile, courant);

  // error file output
  fprintf(stderr, "%4d: t=%.3lfs\n", timestep, Time);

  // prof file output
  if (Time >= settings.outputInterval * double(nfile)) {
    FILE *fp;
    char  filename[256];

    sprintf(filename, "result/prof/output_%04d.prof", nfile);
    fp = fopen(filename, "w");
    fprintf(fp, "%lf\n", Time);
    fprintf(fp, "%d\n", np);
    rep(i, 0, np) {
      if (particles[i].type == ParticleType::Ghost) continue;

      fprintf(fp, "%4d %2d % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 09.3lf % 08.3lf\n",
              i, particles[i].type,
              particles[i].position.x(), particles[i].position.y(), particles[i].position.z(),
              particles[i].velocity[0], particles[i].velocity[1], particles[i].velocity[2],
              P[i], n[i]);
    }
    fclose(fp);

    nfile++;
  }
}

void        cal_gravity() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Fluid) {
      a[i] = settings.gravity;

    } else {
      a[i] << 0.0, 0.0, 0.0;
    }
  }
}

void cal_viscosity() {
  double A = (settings.kinematicViscosity) * (2.0 * settings.dim) / (n0_for_lap * lambda);

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type != ParticleType::Fluid) continue;

    Eigen::Vector3d viscosity_term = Eigen::Vector3d::Zero();

    rep(j_neighbor, 0, num_neighbor[i]) {
      int    j    = neighbor_id[i][j_neighbor];
      double dis2 = neighbor_dis2[i][j_neighbor];

      if (dis2 < re2_for_lap) {
        double dis = sqrt(dis2);
        viscosity_term += (particles[j].velocity - particles[i].velocity) * weight(dis, re_for_lap);
      }
    }

    viscosity_term *= A;
    a[i] += viscosity_term;
  }
}

void        move_particle() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Fluid) {
      particles[i].velocity += a[i] * settings.dt;
      particles[i].position += particles[i].velocity * settings.dt;
    }

    a[i].setZero();
  }
}

void collision() {
  std::vector<Eigen::Vector3d> u_after(np);

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type != ParticleType::Fluid) continue;

    u_after[i] = particles[i].velocity;

    rep(j_neighbor, 0, num_neighbor[i]) {
      int    j    = neighbor_id[i][j_neighbor];
      double dis2 = neighbor_dis2[i][j_neighbor];

      if (dis2 < collision_dis2) {
        double dis     = sqrt(dis2);
        double forceDT = -(particles[j].velocity - particles[i].velocity).dot(particles[j].position - particles[i].position) / dis;  // impulse of collision between particles

        if (forceDT > 0.0) {
          double mi = rho;
          double mj = rho;
          forceDT *= (1.0 + settings.coefficientOfRestitution) * mi * mj / (mi + mj);
          u_after[i] -= (forceDT / mi) * (particles[j].position - particles[i].position) / dis;

#pragma omp critical
          {
            if (j > i) cerr << "WARNING: collision occured between " << i << " and " << j << " particles." << endl;
          }
        }
      }
    }
  }

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type != ParticleType::Fluid) continue;

    particles[i].position += (u_after[i] - particles[i].velocity) * settings.dt;
    particles[i].velocity = u_after[i];
  }
}

void cal_P() {
  cal_n();
  set_boundary_condition();
  set_source_term();
  set_matrix();
  solve_Poisson_eq();
  remove_negative_P();
  set_P_min();
}

void        cal_n() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Ghost) continue;

    n[i] = 0.0;

    rep(j_neighbor, 0, num_neighbor[i]) {
      int    j    = neighbor_id[i][j_neighbor];
      double dis2 = neighbor_dis2[i][j_neighbor];

      if (dis2 < re2_for_n) {
        double dis = sqrt(dis2);
        n[i] += weight(dis, re_for_n);
      }
    }
  }
}

void set_boundary_condition() {
  double n0   = n0_for_n;
  double beta = settings.thresholdForSurfaceDetection;

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Ghost || particles[i].type == ParticleType::DummyWall) {
      boundary_condition[i] = GHOST_OR_DUMMY;

    } else if (n[i] < beta * n0) {
      boundary_condition[i] = SURFACE_PARTICLE;

    } else {
      boundary_condition[i] = INNER_PARTICLE;
    }
  }
}

void set_source_term() {
  double n0    = n0_for_n;
  double gamma = settings.relaxationCoefficientForPressure;

#pragma omp parallel for
  rep(i, 0, np) {
    if (boundary_condition[i] == INNER_PARTICLE) {
      source_term[i] = gamma * (1.0 / (settings.dt * settings.dt)) * ((n[i] - n0) / n0);

    } else {
      source_term[i] = 0.0;
    }
  }
}

void set_matrix() {
  coef_matrix.setZero();

  double n0 = n0_for_lap;
  double A  = 2.0 * settings.dim / (n0 * lambda);
#pragma omp parallel for
  rep(i, 0, np) {
    if (boundary_condition[i] != INNER_PARTICLE) continue;

    rep(j_neighbor, 0, num_neighbor[i]) {
      int j = neighbor_id[i][j_neighbor];
      if (particles[j].type == ParticleType::DummyWall) continue;

      double dis2 = neighbor_dis2[i][j_neighbor];
      if (dis2 < re2_for_lap) {
        double dis     = sqrt(dis2);
        double coef_ij = A * weight(dis, re_for_lap) / rho;

        coef_matrix(i, j) = (-1.0 * coef_ij);
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

void        check_boundary_condition() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (boundary_condition[i] == GHOST_OR_DUMMY) {
      flag_for_checking_boundary_condition[i] = GHOST_OR_DUMMY;

    } else if (boundary_condition[i] == SURFACE_PARTICLE) {
      flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_CONNECTED;

    } else {
      flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  int count;
  while (true) {
    count = 0;

    rep(i, 0, np) {
      if (flag_for_checking_boundary_condition[i] != DIRICHLET_BOUNDARY_IS_CONNECTED) continue;

      rep(j_neighbor, 0, num_neighbor[i]) {
        int j = neighbor_id[i][j_neighbor];
        if (flag_for_checking_boundary_condition[j] != DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) continue;

        double dis2 = neighbor_dis2[i][j_neighbor];
        if (dis2 < re2_for_lap) flag_for_checking_boundary_condition[j] = DIRICHLET_BOUNDARY_IS_CONNECTED;
      }

      flag_for_checking_boundary_condition[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
      count++;
    }

    if (count == 0) break;
  }

#pragma omp parallel for
  rep(i, 0, np) {
    if (flag_for_checking_boundary_condition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
#pragma omp critical
      {
        cerr << "WARNING: There is no dirichlet boundary condition for particle " << i << endl;
      }
    }
  }
}

void        increase_diagonal_term() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (flag_for_checking_boundary_condition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
      coef_matrix(i, i) *= 2.0;
    }
  }
}

void solve_Poisson_eq() {
  Eigen::SparseMatrix<double>                  A = coef_matrix.sparseView();
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  P = solver.solveWithGuess(source_term, P);
}

void        remove_negative_P() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (P[i] < 0.0) P[i] = 0.0;
  }
}

void        set_P_min() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (boundary_condition[i] == GHOST_OR_DUMMY) continue;

    P_min[i] = P[i];

    rep(j_neighbor, 0, num_neighbor[i]) {
      int j = neighbor_id[i][j_neighbor];
      if (particles[j].type == ParticleType::DummyWall) continue;

      double dis2 = neighbor_dis2[i][j_neighbor];
      if (dis2 < re2_for_grad) {
        if (P_min[i] > P[j]) P_min[i] = P[j];
      }
    }
  }
}

void cal_P_grad() {
  double A = settings.dim / n0_for_grad;

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type != ParticleType::Fluid) continue;

    Eigen::Vector3d grad = Eigen::Vector3d::Zero();

    rep(j_neighbor, 0, num_neighbor[i]) {
      int j = neighbor_id[i][j_neighbor];
      if (particles[j].type == ParticleType::DummyWall) continue;

      double dis2 = neighbor_dis2[i][j_neighbor];
      if (dis2 < re2_for_grad) {
        double dis = sqrt(dis2);
        grad += (particles[j].position - particles[i].position) * (P[j] - P_min[i]) * weight(dis, re_for_grad) / dis2;
      }
    }

    grad *= A;
    a[i] = -1.0 * grad / rho;
  }
}

void        move_particle_using_P_grad() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Fluid) {
      particles[i].velocity += a[i] * settings.dt;
      particles[i].position += a[i] * settings.dt * settings.dt;
    }

    a[i].Zero();
  }
}

void cal_courant() {
  courant = 0.0;

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type != ParticleType::Fluid) continue;

    double courant_i = (particles[i].velocity.norm() * settings.dt) / settings.particleDistance;
    if (courant_i > courant)
#pragma omp critical
    {
      courant = courant_i;
    }
  }

  if (courant > settings.cflCondition) {
    cerr << "ERROR: Courant number is larger than CFL condition. Courant = " << courant << endl;
    error_flag = ON;
  }
}

void        store_particle() {
#pragma omp parallel for
  rep(i, 0, num_bucket) {
    bucket_first[i] = -1;
    bucket_last[i]  = -1;
  }

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Ghost) continue;
    bucket_next[i] = -1;
  }

#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Ghost) continue;

    int ix      = (int)((particles[i].position.x() - x_min) / bucket_length) + 1;
    int iy      = (int)((particles[i].position.y() - y_min) / bucket_length) + 1;
    int iz      = (int)((particles[i].position.z() - z_min) / bucket_length) + 1;
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

void        set_neighbor() {
#pragma omp parallel for
  rep(i, 0, np) {
    if (particles[i].type == ParticleType::Ghost) continue;

    num_neighbor[i] = 0;

    int ix = int((particles[i].position.x() - x_min) / bucket_length) + 1;
    int iy = int((particles[i].position.y() - y_min) / bucket_length) + 1;
    int iz = int((particles[i].position.z() - z_min) / bucket_length) + 1;

    for (int jx = ix - 1; jx <= ix + 1; jx++) {
      for (int jy = iy - 1; jy <= iy + 1; jy++) {
        for (int jz = iz - 1; jz <= iz + 1; jz++) {
          int jbucket = jz * num_bucket_xy + jy * num_bucket_x + jx;
          int j       = bucket_first[jbucket];

          while (j != -1) {
            double dis2 = cal_dis2(i, j);
            if (j != i && dis2 < re2_max) {
              if (num_neighbor[i] + 1 == NEIGHBOR_ARRAY_SIZE) {
                cerr << "ERROR: Neighbor size of particle " << i << " is bigger than expected." << endl;
                cout << "ERROR: Neighbor size of particle " << i << " is bigger than expected." << endl;
                exit(1);
              }

              neighbor_id[i][num_neighbor[i]]   = j;
              neighbor_dis2[i][num_neighbor[i]] = dis2;
              num_neighbor[i]++;
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

  if (dis < re) w = (re / dis) - 1.0;

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