#include "simulation.hpp"

#include "bucket.hpp"
#include "csv.hpp"
#include "mps.hpp"
#include "particle.hpp"
#include "settings.hpp"

#include <Eigen/Dense>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

#define rep(i, a, b) for (int i = a; i < b; i++)

void Simulation::run() {
    startSimulation();

    std::vector<Particle> particles;

    read_data(particles);
    mps = MPS(settings, particles);

    timestep = 0;

    write_data(0.0, system_clock::now());

    startTime = system_clock::now();
    while (time <= settings.finishTime) {
        system_clock::time_point timestepStartTime = system_clock::now();

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

    logFile.open("result/result.csv");
    if (!logFile.is_open()) {
        cerr << "ERROR: Could not open the log file: " << resultFileNum << std::endl;
        exit(-1);
    }
    auto logFileWriter = csv::make_csv_writer(logFile);
    logFileWriter << std::vector<std::string>{
        "Timestep",
        "Dt (s)",
        "Current Time (s)",
        "Finish Time (s)",
        "Elapsed Time (h:m:s)",
        "Elapsed Time (s)",
        "Remaining Time (h:m:s)",
        "Average Time per Timestep (s)",
        "Time for This Timestep (s)",
        "Number of Output Files",
        "Courant Number"
    };
}

void Simulation::endSimulation() {
    cout << endl
         << std::format(
                "Total Simulation Time = {:%Hh %Mm %Ss}",
                duration_cast<seconds>(system_clock::now() - startTime)
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
    const double& courantNumber, const system_clock::time_point& timestepStartTime
) {
    seconds elapsed = duration_cast<seconds>(system_clock::now() - startTime);

    double average = 0;
    if (timestep != 0) {
        average = duration_cast<milliseconds>(system_clock::now() - startTime).count() /
                  (double) (1000 * timestep);
    }

    seconds remain{int(((settings.finishTime - time) / time) * average * timestep)};

    // last
    milliseconds last =
        duration_cast<milliseconds>(system_clock::now() - timestepStartTime);

    std::string formattedTime          = std::format("{:.3f}", time);
    std::string formattedElapsed       = std::format("{:%Hh %Mm %Ss}", elapsed);
    std::string formattedAverage       = std::format("{:.3f}", average);
    std::string formattedRemain        = std::format("{:%Hh %Mm %Ss}", remain);
    std::string formattedLast          = std::format("{:%S}", last);
    std::string formattedCourantNumber = std::format("{:.2f}", courantNumber);

    cout << std::format(
                "{}: dt={}s   t={}s   fin={}s   elapsed={}   remain={}   "
                "ave={}s/step   last={}s/step   out={}files   Courant={}",
                timestep,
                settings.dt,
                formattedTime,
                settings.finishTime,
                formattedElapsed,
                formattedRemain,
                formattedAverage,
                formattedLast,
                resultFileNum,
                formattedCourantNumber
            )
         << endl;

    auto logFileWriter = csv::make_csv_writer(logFile);
    logFileWriter << std::make_tuple(
        timestep,
        settings.dt,
        formattedTime,
        settings.finishTime,
        formattedElapsed,
        elapsed.count(),
        formattedRemain,
        formattedAverage,
        formattedLast,
        resultFileNum,
        formattedCourantNumber
    );

    if (time >= settings.outputInterval * double(resultFileNum)) {
        std::string filename =
            std::format("result/prof/output_{:04d}.prof", resultFileNum);

        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            throw std::runtime_error(
                std::format("Could not open result file: {}", filename)
            );
        }

        outFile << time << endl;
        outFile << mps.particles.size() << endl;
        for (auto& pi : mps.particles) {
            if (pi.type == ParticleType::Ghost) {
                continue;
            }

            outFile << std::format(
                           "{:4d} {:2d} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} "
                           "{:9.3f} "
                           "{:8.3f}",
                           pi.id,
                           static_cast<int>(pi.type),
                           pi.position.x(),
                           pi.position.y(),
                           pi.position.z(),
                           pi.velocity.x(),
                           pi.velocity.y(),
                           pi.velocity.z(),
                           pi.pressure,
                           pi.numberDensity
                       )
                    << endl;
        }

        outFile.close();

        resultFileNum++;
    }
}
