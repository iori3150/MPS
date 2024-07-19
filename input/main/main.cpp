#include <bits/stdc++.h>
using namespace std;
#define rep(i, a, b) for (int i = a; i < b; i++)

#define PARTICLE_DISTANCE 0.05

// type
#define GHOST -1
#define FLUID 0
#define WALL 2
#define DUMMY_WALL 3

#define EPS (0.01 * PARTICLE_DISTANCE)

int            np;
vector<int>    type;
vector<double> x;

double l0 = PARTICLE_DISTANCE;

// common
int itype;

// set_line
vector<double> x_range(2), y_range(2), z_range(2);

// cal_min_and_max
double x_min, x_max, y_min, y_max, z_min, z_max;

void set_particle();
void cal_min_and_max();
void check_closeness();
void write_data();

int main() {
  np = 0;

  double left   = -1.0;
  double right  = 1.0;
  double bottom = 0.0;
  double top    = 1.0;

  // *****************
  //
  // ***** fluid *****
  //
  // *****************
  itype   = FLUID;
  x_range = {left + l0 / 2.0, left + (right - left) / 4.0 - l0 / 2.0};
  y_range = {bottom + l0 / 2.0, top * 0.8 - l0 / 2.0};
  z_range = {0.0, 0.0};
  set_particle();

  // *********************
  //
  // ***** left wall *****
  //
  // *********************
  itype   = WALL;
  x_range = {left - l0 / 2.0, left - l0 / 2.0};
  y_range = {bottom + l0 / 2.0, top - l0 / 2.0};
  z_range = {0.0, 0.0};
  set_particle();
  itype   = DUMMY_WALL;
  x_range = {left - l0 / 2.0 - 3.0 * l0, left - l0 / 2.0 - l0};
  // no change in y_range and z_range
  set_particle();

  // **********************
  //
  // ***** right wall *****
  //
  // **********************
  itype   = WALL;
  x_range = {right + l0 / 2.0, right + l0 / 2.0};
  y_range = {bottom + l0 / 2.0, top - l0 / 2.0};
  z_range = {0.0, 0.0};
  set_particle();
  itype   = DUMMY_WALL;
  x_range = {right + l0 / 2.0 + l0, right + l0 / 2.0 + 3.0 * l0};
  // no change in y_range and z_range
  set_particle();

  // ***********************
  //
  // ***** bottom wall *****
  //
  // ***********************
  itype   = WALL;
  x_range = {left - l0 / 2.0, right + l0 / 2.0};
  y_range = {bottom - l0 / 2.0, bottom - l0 / 2.0};
  z_range = {0.0, 0.0};
  set_particle();
  itype   = DUMMY_WALL;
  y_range = {bottom - l0 / 2.0 - 3.0 * l0, bottom - l0 / 2.0 - l0};
  // no change in x_range and z_range
  set_particle();

  // *************************
  //
  // ***** bottom corner *****
  //
  // *************************
  itype   = DUMMY_WALL;
  x_range = {left - l0 / 2.0 - 3.0 * l0, left - l0 / 2.0 - l0};
  y_range = {bottom - l0 / 2.0 - 3.0 * l0, bottom - l0 / 2.0};
  z_range = {0.0, 0.0};
  set_particle();
  x_range = {right + l0 / 2.0 + l0, right + l0 / 2.0 + 3 * l0};
  set_particle();

  cal_min_and_max();

  check_closeness();

  write_data();

  return 0;
}

void set_particle() {
  int ix, iy, iz;

  ix = 0;
  while (x_range[0] + double(ix) * l0 <= x_range[1] + EPS) {
    iy = 0;
    while (y_range[0] + double(iy) * l0 <= y_range[1] + EPS) {
      iz = 0;
      while (z_range[0] + double(iz) * l0 <= z_range[1] + EPS) {
        type.push_back(itype);
        x.push_back(x_range[0] + double(ix) * l0);
        x.push_back(y_range[0] + double(iy) * l0);
        x.push_back(z_range[0] + double(iz) * l0);
        np++;

        iz++;
      }
      iy++;
    }
    ix++;
  }
}

void cal_min_and_max() {
  rep(i, 0, np) {
    if (type[i] == GHOST) continue;
    x_min = x_max = x[i * 3 + 0];
    y_min = y_max = x[i * 3 + 1];
    z_min = z_max = x[i * 3 + 2];
    break;
  }

  rep(i, 0, np) {
    if (type[i] == GHOST) continue;

    if (x_min > x[i * 3 + 0]) x_min = x[i * 3 + 0];
    if (x_max < x[i * 3 + 0]) x_max = x[i * 3 + 0];
    if (y_min > x[i * 3 + 1]) y_min = x[i * 3 + 1];
    if (y_max < x[i * 3 + 1]) y_max = x[i * 3 + 1];
    if (z_min > x[i * 3 + 2]) z_min = x[i * 3 + 2];
    if (z_max < x[i * 3 + 2]) z_max = x[i * 3 + 2];
  }
}

void check_closeness() {
  double x_ij, y_ij, z_ij;
  double dis;
  rep(i, 0, np) {
    if (type[i] == GHOST) continue;
    rep(j, 0, np) {
      if (j == i || type[j] == GHOST) continue;

      x_ij = x[j * 3 + 0] - x[i * 3 + 0];
      y_ij = x[j * 3 + 1] - x[i * 3 + 1];
      z_ij = x[j * 3 + 2] - x[i * 3 + 2];
      dis  = sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij);
      if (dis < 0.1 * l0 && i < j) {
        cout << "WARNING: There are too close particles. Check error.log" << endl;

        cerr << "WARNING: Particle " << i << " and " << j << " are too close." << endl;
        cerr << "         type[" << i << "]= " << type[i] << " x[" << i << "] = (" << x[i * 3] << ", " << x[i * 3 + 1] << ", " << x[i * 3 + 2] << ")" << endl;
        cerr << "         type[" << j << "]= " << type[j] << " x[" << j << "] = (" << x[j * 3] << ", " << x[j * 3 + 1] << ", " << x[j * 3 + 2] << ")" << endl;
      }
    }
  }
}

void write_data() {
  FILE *fp;
  char  filename[256];

  sprintf(filename, "result/input.prof");
  fp = fopen(filename, "w");
  fprintf(fp, "%lf\n", 0.0);
  fprintf(fp, "%d\n", np);
  rep(i, 0, np)
      fprintf(fp, "%4d %2d % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf % 08.3lf\n", i, type[i], x[i * 3], x[i * 3 + 1], x[i * 3 + 2], 0.0, 0.0, 0.0, 0.0, 0.0);
  fclose(fp);

  sprintf(filename, "result/input.data");
  fp = fopen(filename, "w");
  fprintf(fp, "x_min %lf\n", x_min);
  fprintf(fp, "x_max %lf\n", x_max);
  fprintf(fp, "y_min %lf\n", y_min);
  fprintf(fp, "y_max %lf\n", y_max);
  fprintf(fp, "z_min %lf\n", z_min);
  fprintf(fp, "z_max %lf\n", z_max);
  fclose(fp);
}
