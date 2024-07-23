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

    saver.save(mps.particles, time, resultFileNum);

    simulationStartTime = system_clock::now();
    while (time <= settings.finishTime) {
        auto timestepStartTime = system_clock::now();

        mps.stepForward();
        timestep++;
        time += settings.dt;

        auto timestepEndTime = system_clock::now();

        if (time >= settings.outputInterval * double(resultFileNum)) {
            saver.save(mps.particles, time, resultFileNum);
            resultFileNum++;
        }
        timeStepReport(timestepStartTime, timestepEndTime, mps.getCourantNumber());
    }
    simulationEndTime = system_clock::now();

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

    saver = Saver("result");
}

void Simulation::endSimulation() {
    auto totalSimulationTime =
        duration_cast<seconds>(simulationEndTime - simulationStartTime);

    cout << endl
         << std::format("Total Simulation Time = {:%Hh %Mm %Ss}", totalSimulationTime)
         << endl;

    cout << endl << "*** END SIMULATION ***" << endl << endl;

    logFile.close();
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

void Simulation::timeStepReport(
    const system_clock::time_point& timeStepStartTime,
    const system_clock::time_point& timeStepEndTime,
    const double& courantNumber
) {
    auto elapsed = duration_cast<seconds>(timeStepEndTime - simulationStartTime);

    double average = 0;
    if (timestep != 0) {
        average =
            duration_cast<milliseconds>(timeStepEndTime - simulationStartTime).count() /
            (double) (1000 * timestep);
    }

    auto remain{int(((settings.finishTime - time) / time) * average * timestep)};

    auto last = duration_cast<milliseconds>(timeStepEndTime - timeStepStartTime);

    auto formattedTime          = std::format("{:.3f}", time);
    auto formattedElapsed       = std::format("{:%Hh %Mm %Ss}", elapsed);
    auto formattedAverage       = std::format("{:.3f}", average);
    auto formattedRemain        = std::format("{:%Hh %Mm %Ss}", remain);
    auto formattedLast          = std::format("{:%S}", last);
    auto formattedCourantNumber = std::format("{:.2f}", courantNumber);

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
}
