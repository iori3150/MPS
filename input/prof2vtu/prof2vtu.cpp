#include <bits/stdc++.h>
using namespace std;
#define rep(i, a, b) for (int i = a; i < b; i++)

#define IN_FILE "result/input.prof"
#define OUT_FILE "result/input.vtu"

int            np;  // number of particles
vector<int>    type;
vector<double> x, u, P, n;

void read_data();
void write_data();

int main(void) {
  cout << "*** START CONVERTING ***" << endl;

  read_data();

  write_data();

  cout << "*** END CONVERTING *** " << endl;

  return 0;
}

void read_data() {
  ifstream file;
  file.open(IN_FILE);

  if (!file) {
    cout << "ERROR: There is no file named " << IN_FILE << endl;
    exit(1);
  }

  double Time;
  file >> Time;
  file >> np;

  type.resize(np);
  x.resize(np * 3);
  u.resize(np * 3);
  P.resize(np);
  n.resize(np);

  int id;
  rep(i, 0, np) {
    file >> id >> type[i];
    file >> x[i * 3] >> x[i * 3 + 1] >> x[i * 3 + 2];
    file >> u[i * 3] >> u[i * 3 + 1] >> u[i * 3 + 2];
    file >> P[i] >> n[i];
  }

  file.close();
}

void write_data() {
  char filename[256];
  sprintf(filename, OUT_FILE);

  cout << filename << endl;

  FILE *fp;
  fp = fopen(filename, "w");

  fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", np, np);

  // Position
  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%lf %lf %lf\n", x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<PointData>\n");

  // ParticleType
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%d\n", type[i]);
  fprintf(fp, "</DataArray>\n");

  // Velocity
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity_abs' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", sqrt(u[i * 3] * u[i * 3] + u[i * 3 + 1] * u[i * 3 + 1] + u[i * 3 + 2] * u[i * 3 + 2]));
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity_x' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity_y' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3 + 1]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity_z' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3 + 2]);
  fprintf(fp, "</DataArray>\n");

  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%d\n", i);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%d\n", i + 1);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "1\n");
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);
}