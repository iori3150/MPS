#include "simulation.hpp"

#include "bucket.hpp"
#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <Eigen/Dense>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;
namespace chrono = std::chrono;
namespace fs     = std::filesystem;

#define rep(i, a, b) for (int i = a; i < b; i++)

void Simulation::run() {
    startSimulation();

    std::vector<Particle> particles;

    read_data(particles);
    mps = MPS(settings, particles);

    timestep = 0;

    write_data(0.0, chrono::system_clock::now());

    while (time <= settings.finishTime) {
        chrono::system_clock::time_point timestepStartTime = chrono::system_clock::now();

        mps.stepForward();
        double courantNumber = mps.calcCourantNumber();

        timestep++;
        time += settings.dt;
        write_data(courantNumber, timestepStartTime);
    }

    endSimulation();
}

void Simulation::startSimulation() {
    cout << endl << "*** START SIMULATION ***" << endl;

    startTime = chrono::system_clock::now();

    logFile.open("result/result.log");
    if (!logFile.is_open()) {
        cerr << "ERROR: Could not open the log file: " << resultFileNum << std::endl;
        exit(-1);
    }
}

void Simulation::endSimulation() {
    cout << endl
         << std::format(
                "Total Simulation Time = {:%Hh %Mm %Ss}",
                chrono::duration_cast<chrono::seconds>(
                    chrono::system_clock::now() - startTime
                )
            )
         << endl;

    logFile.close();

    cout << endl << "*** END SIMULATION ***" << endl << endl;
}

void Simulation::read_data(std::vector<Particle>& particles) {
    std::ifstream file;

    file.open(settings.inputProfPath);
    if (!file) {
        cout << "ERROR: There is no file named " << settings.inputProfPath << endl;
        exit(1);
    }

    int numberOfParticles;
    file >> time;
    file >> numberOfParticles;

    int id;
    rep(i, 0, numberOfParticles) {
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
    settings.domain = Domain(xMin, xMax, yMin, yMax, zMin, zMax);

    file.close();
}

void Simulation::write_data(
    const double& courantNumber, const chrono::system_clock::time_point& timestepStartTime
) {
    // elapsed
    char elapsed[256];
    sprintf(
        elapsed,
        "elapsed=%s",
        std::format("{:%Hh %Mm %Ss}", chrono::system_clock::now() - startTime).c_str()
    );

    // ave [s]/[timestep]
    double ave = 0;
    if (timestep != 0)
        ave = chrono::duration_cast<chrono::seconds>(
                  chrono::system_clock::now() - startTime
              )
                  .count() /
              timestep;

    // remain
    char remain[256];
    int remainingTimeInt = ((settings.finishTime - time) / time) * ave * timestep;
    chrono::seconds remainingTime{remainingTimeInt};
    if (timestep == 0)
        sprintf(remain, "remain=---");
    else
        sprintf(
            remain,
            "remain=%s",
            std::format("{:%Hh %Mm %Ss}", remainingTime).c_str()
        );

    // last
    chrono::milliseconds last = chrono::duration_cast<chrono::milliseconds>(
        chrono::system_clock::now() - timestepStartTime
    );

    // terminal output
    printf(
        "%d: settings.dt=%.gs   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs/step   "
        "last=%.3lfs/step   out=%dfiles   Courant=%.2lf\n",
        timestep,
        settings.dt,
        time,
        settings.finishTime,
        elapsed,
        remain,
        ave,
        last,
        resultFileNum,
        courantNumber
    );

    logFile << std::format(
                   "{}: dt={}s   t={:.3f}s   fin={:.1f}s   {}   {}   ave={:.3f}s/step "
                   "last={:%S}s/step   out={}files   Courant={:.2f}",
                   timestep,
                   settings.dt,
                   time,
                   settings.finishTime,
                   elapsed,
                   remain,
                   ave,
                   last,
                   resultFileNum,
                   courantNumber
               )
            << endl;

    // error file output
    fprintf(stderr, "%4d: t=%.3lfs\n", timestep, time);

    // prof file output
    if (time >= settings.outputInterval * double(resultFileNum)) {
        FILE* fp;
        char filename[256];

        sprintf(filename, "result/prof/output_%04d.prof", resultFileNum);
        fp = fopen(filename, "w");
        fprintf(fp, "%lf\n", time);
        fprintf(fp, "%d\n", mps.particles.size());
        rep(i, 0, mps.particles.size()) {
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

        resultFileNum++;
    }
}

std::tuple<int, int, int> Simulation::getTimeDuration(
    const chrono::system_clock::time_point& start,
    const chrono::system_clock::time_point& end
) {
    chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
    chrono::hours hours      = chrono::duration_cast<chrono::hours>(duration);
    chrono::minutes minutes  = chrono::duration_cast<chrono::minutes>(duration - hours);
    chrono::seconds seconds  = duration - hours - minutes;
    return std::forward_as_tuple(hours.count(), minutes.count(), seconds.count());
}
