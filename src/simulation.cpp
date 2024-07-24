#include "simulation.hpp"

#include "bucket.hpp"
#include "csv.hpp"
#include "exporter.hpp"
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
namespace fs = std::filesystem;

void Simulation::run() {
    startSimulation();

    std::vector<Particle> particles;

    read_data(particles);
    mps = MPS(settings, particles);

    exportParticles(mps.particles);

    simulationStartTime = system_clock::now();
    while (time <= settings.finishTime) {
        auto timeStepStartTime = system_clock::now();

        time += settings.dt;
        mps.stepForward(isTimeToExport());
        if (isTimeToExport()) {
            exportParticles(mps.particles);
        }

        auto timeStepEndTime = system_clock::now();
        timeStepReport(timeStepStartTime, timeStepEndTime, mps.getCourantNumber());

        timeStep++;
    }
    simulationEndTime = system_clock::now();

    endSimulation();
}

void Simulation::startSimulation() {
    cout << endl << "*** START SIMULATION ***" << endl;

    logFile.open("result/log.csv");
    if (!logFile.is_open()) {
        cout << "ERROR: Could not open the log file: "
             << "result/log.csv" << std::endl;
        exit(-1);
    }
    auto logFileWriter = csv::make_csv_writer(logFile);
    logFileWriter << std::vector<std::string>{
        "Time Step",
        "Dt (s)",
        "Current Time (s)",
        "Finish Time (s)",
        "Elapsed Time (h:m:s)",
        "Elapsed Time (s)",
        "Remaining Time (h:m:s)",
        "Average Time per Time Step (s)",
        "Time for This Time Step (s)",
        "Number of Output Files",
        "Courant Number"
    };
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
    int particleDataHeaderRow = 3;

    // Set up CSV format for meta data
    csv::CSVFormat metaDataFormat;
    metaDataFormat.no_header();
    csv::CSVReader metaDataReader(settings.inputProfPath, metaDataFormat);

    // Read meta data
    for (auto& row : metaDataReader) {
        // Get the time from the first row
        if (metaDataReader.n_rows() + 1 == 1)
            time = row[1].get<double>();

        // Stop reading after the row before the particle data header
        if (metaDataReader.n_rows() + 1 == particleDataHeaderRow - 1)
            break;
    }

    // Set up CSV format for particle data
    csv::CSVFormat particleDataFormat;
    particleDataFormat.header_row(particleDataHeaderRow - 1);
    csv::CSVReader particleDataReader(settings.inputProfPath, particleDataFormat);

    // Read particle data and create Particle objects
    for (auto& row : particleDataReader) {
        int id      = row["ID"].get<int>();
        double type = row["Type"].get<int>();
        double x    = row["Position.x (m)"].get<double>();
        double y    = row["Position.y (m)"].get<double>();
        double z    = row["Position.z (m)"].get<double>();
        double u    = row["Velocity.x (m/s)"].get<double>();
        double v    = row["Velocity.y (m/s)"].get<double>();
        double w    = row["Velocity.z (m/s)"].get<double>();

        particles.push_back(Particle(
            id,
            static_cast<ParticleType>(type),
            Eigen::Vector3d(x, y, z),
            Eigen::Vector3d(u, v, w),
            settings.density
        ));
    }
    // std::ifstream file;
    // file.open(settings.inputDataPath);
    // if (!file) {
    //     cout << "ERROR: There is no file named " << settings.inputDataPath << endl;
    //     exit(1);
    // }

    // std::string dummy_string;
    // double xMin, xMax, yMin, yMax, zMin, zMax;
    // file >> dummy_string >> xMin;
    // file >> dummy_string >> xMax;
    // file >> dummy_string >> yMin;
    // file >> dummy_string >> yMax;
    // file >> dummy_string >> zMin;
    // file >> dummy_string >> zMax;
    // settings.domain = Domain(xMin, xMax, yMin, yMax, zMin, zMax);

    // file.close();
    settings.domain = Domain(-0.1, 1.1, -0.1, 0.8, -1, 1);
}

void Simulation::timeStepReport(
    const system_clock::time_point& timeStepStartTime,
    const system_clock::time_point& timeStepEndTime,
    const double& courantNumber
) {
    auto elapsed = duration_cast<seconds>(timeStepEndTime - simulationStartTime);

    double average = 0;
    if (timeStep != 0) {
        average =
            duration_cast<milliseconds>(timeStepEndTime - simulationStartTime).count() /
            (double) (1000 * timeStep);
    }

    seconds remain{int(((settings.finishTime - time) / time) * average * timeStep)};

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
                timeStep,
                settings.dt,
                formattedTime,
                settings.finishTime,
                formattedElapsed,
                formattedRemain,
                formattedAverage,
                formattedLast,
                outFileNum,
                formattedCourantNumber
            )
         << endl;

    auto logFileWriter = csv::make_csv_writer(logFile);
    logFileWriter << std::make_tuple(
        timeStep,
        settings.dt,
        formattedTime,
        settings.finishTime,
        formattedElapsed,
        elapsed.count(),
        formattedRemain,
        formattedAverage,
        formattedLast,
        outFileNum,
        formattedCourantNumber
    );
}

bool Simulation::isTimeToExport() {
    return time >= settings.outputInterval * double(outFileNum);
}

void Simulation::exportParticles(const std::vector<Particle>& particles) {
    Exporter exporter;
    exporter.toCsv(
        fs::path(std::format("result/csv/output_{:04}.csv", outFileNum)),
        particles,
        time
    );
    exporter.toVtu(
        fs::path(std::format("result/vtu/output_{:04}.vtu", outFileNum)),
        particles,
        time
    );

    outFileNum++;
}
