#include "simulation.hpp"

#include "bucket.hpp"
#include "csv.hpp"
#include "exporter.hpp"
#include "mps.hpp"
#include "particle.hpp"

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
namespace fs     = std::filesystem;
namespace chrono = std::chrono;

void Simulation::run() {
    startSimulation();

    mps  = MPS();
    time = mps.initialize();
    exportParticles(mps.particles);

    simulationStartTime = system_clock::now();
    while (time <= mps.settings.finishTime) {
        auto timeStepStartTime = system_clock::now();

        time += mps.settings.dt;
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
} // namespace std::chrono

void Simulation::startSimulation() {
    cout << endl << "*** START SIMULATION ***" << endl;

    createResultDirectory();

    logFile.open(resultDirectory / "log.csv");
    if (!logFile.is_open()) {
        cout << "ERROR: Could not open the log file: " << resultDirectory / "log.csv"
             << std::endl;
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

void Simulation::createResultDirectory() {
    fs::path parentDirectory = "./results";
    fs::create_directory(parentDirectory);

    const auto currentTime =
        chrono::zoned_time{chrono::current_zone(), system_clock::now()};
    resultDirectory = parentDirectory / format("{:%F_%H-%M}", currentTime);
    int count       = 1;
    while (fs::exists(resultDirectory)) {
        resultDirectory = parentDirectory / format("{:%F_%H-%M}({})", currentTime, count);
        count++;
    }
    fs::create_directory(resultDirectory);
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

    seconds remain{int(((mps.settings.finishTime - time) / time) * average * timeStep)};

    auto last = duration_cast<milliseconds>(timeStepEndTime - timeStepStartTime);

    auto formattedTime          = std::format("{:.3f}", time);
    auto formattedElapsed       = std::format("{:%H:%M:%S}", elapsed);
    auto formattedAverage       = std::format("{:.3f}", average);
    auto formattedRemain        = std::format("{:%H:%M:%S}", remain);
    auto formattedLast          = std::format("{:%S}", last);
    auto formattedCourantNumber = std::format("{:.2f}", courantNumber);

    cout << std::format(
                "{}: dt={}s   t={}s   fin={}s   elapsed={}   remain={}   "
                "ave={}s/step   last={}s/step   out={}files   Courant={}",
                timeStep,
                mps.settings.dt,
                formattedTime,
                mps.settings.finishTime,
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
        mps.settings.dt,
        formattedTime,
        mps.settings.finishTime,
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
    return time >= mps.settings.outputInterval * double(outFileNum);
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
