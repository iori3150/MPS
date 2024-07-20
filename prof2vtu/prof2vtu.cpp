#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
#define rep(i, a, b) for (int i = a; i < b; i++)

#define IN_DIR "result/prof/"
#define OUT_DIR "result/vtu/"

int            np;  // number of particles
vector<int>    type;
vector<double> x, u, P, n;

int nfile;

void read_data(string path);
void write_data();

int main(void) {
  cout << "*** START CONVERTING ***" << endl;

  nfile = 0;
  for (const auto& file : filesystem::directory_iterator(IN_DIR)) {
    read_data(file.path().string());
    write_data();

    nfile++;
  }

  cout << "*** END CONVERTING *** " << endl;

  return 0;
}

void read_data(string path) {
  ifstream file;
  file.open(path);

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
  sprintf(filename, "%soutput_%05d.vtu", OUT_DIR, nfile);

  cout << filename << endl;

  FILE* fp;
  fp = fopen(filename, "w");

  fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", np, np);

  // position
  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%lf %lf %lf\n", x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<PointData>\n");

  // particle type
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Int32' Name='Particle Type' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%d\n", type[i]);
  fprintf(fp, "</DataArray>\n");

  // velocity
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", sqrt(u[i * 3] * u[i * 3] + u[i * 3 + 1] * u[i * 3 + 1] + u[i * 3 + 2] * u[i * 3 + 2]));
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity X' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity Y' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3 + 1]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity Z' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", u[i * 3 + 2]);
  fprintf(fp, "</DataArray>\n");

  // pressure
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", P[i]);
  fprintf(fp, "</DataArray>\n");

  // number density
  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Number Density' format='ascii'>\n");
  rep(i, 0, np) fprintf(fp, "%f\n", n[i]);
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