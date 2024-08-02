#include "simulation.hpp"

#include "bucket.hpp"
#include "csv.hpp"
#include "exporter.hpp"
#include "mps.hpp"
#include "particle.hpp"
#include "utilities.hpp"

#include <Eigen/Dense>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
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

void Simulation::run(const fs::path& inputYamlPath) {
    createResultDirectory(inputYamlPath);
    prepareLogFile();

    mps  = MPS(inputYamlPath);
    time = mps.loadInitialState();
    exportParticles(mps.particles);

    spdlog::info("START SIMULATION");
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
}

void Simulation::createResultDirectory(const fs::path& inputYamlPath) {
    fs::path parentDirectory = "./results";
    fs::create_directory(parentDirectory);

    const auto currentTime =
        chrono::zoned_time{chrono::current_zone(), system_clock::now()};
    std::string timestamp = format("{:%F_%H-%M}", currentTime);

    resultDirectory = parentDirectory / timestamp;

    // If the result folder with the same timestamp already exists (i.e. previous
    // simulations have run within one minute), create a new folder name by appending a
    // number. For example, if "2024-08-01_10-02" already exists, create
    // "2024-08-01_10-02(1)", "---(2)", "---(3)", ... until the name is unique.
    int count = 1;
    while (fs::exists(resultDirectory)) {
        resultDirectory = parentDirectory / format("{}({})", timestamp, count);
        count++;
    }

    fs::create_directory(resultDirectory);
    fs::create_directory(resultDirectory / "csv");
    fs::create_directory(resultDirectory / "vtu");

    fs::copy(inputYamlPath, resultDirectory / inputYamlPath.filename());
}

void Simulation::prepareLogFile() {
    // Create console sink and file sink
    auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto file_sink   = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
        fs::path(resultDirectory / "execution.log").string(),
        true
    );

    // Create logger by ombining console sink and file sink
    std::vector<spdlog::sink_ptr> sinks{stdout_sink, file_sink};
    auto logger =
        std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

    // Register logger
    spdlog::set_default_logger(logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");

    fs::path logFilePath = resultDirectory / "log.csv";
    logFile.open(logFilePath);
    if (!logFile.is_open()) {
        exitWithError("Could not open the log file: " + logFilePath.string());
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

    spdlog::info("END SIMULATION");
    spdlog::info(std::format("Total Simulation Time = {:%T}", totalSimulationTime));

    logFile.close();
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
    auto formattedElapsed       = std::format("{:%T}", elapsed);
    auto formattedAverage       = std::format("{:.3f}", average);
    auto formattedRemain        = std::format("{:%T}", remain);
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
        resultDirectory / std::format("csv/output_{:04}.csv", outFileNum),
        particles,
        time
    );
    exporter.toVtu(
        resultDirectory / std::format("vtu/output_{:04}.vtu", outFileNum),
        particles,
        time
    );

    outFileNum++;
}
