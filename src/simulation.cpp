#include "simulation.hpp"

#include "bucket.hpp"
#include "exporter.hpp"
#include "mps.hpp"
#include "particle.hpp"

#include <Eigen/Dense>
#include <cmath>   // for fmod
#include <cstdlib> // for std::exit
#include <csv.hpp>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>

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
    while (time + mps.settings.dt <= mps.settings.finishTime) {
        auto timeStepStartTime   = system_clock::now();
        int initialDebugLogCount = mps.debugLogCount;

        time += mps.settings.dt;

        double courantNumber = mps.stepForward(isTimeToExport());

        if (isTimeToExport()) {
            exportParticles(mps.particles);
        }

        auto timeStepEndTime = system_clock::now();
        timeStepReport(
            timeStepStartTime,
            timeStepEndTime,
            courantNumber,
            mps.debugLogCount > initialDebugLogCount
        );

        timeStep++;
    }
    simulationEndTime = system_clock::now();

    endSimulation();
}

void Simulation::createResultDirectory(const fs::path& inputYamlPath) {
    resultDirectory = inputYamlPath.parent_path() / "result";
    if (fs::exists(resultDirectory)) {
        fs::remove_all(resultDirectory);
    }
    fs::create_directory(resultDirectory);

    fs::create_directory(resultDirectory / "csv");
    fs::create_directory(resultDirectory / "vtu");

    fs::copy(inputYamlPath, resultDirectory / inputYamlPath.filename());
}

void Simulation::prepareLogFile() {
    // Create console sink
    auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    stdout_sink->set_level(spdlog::level::info);

    // Create file sink
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
        fs::path(resultDirectory / "execution.log").string(),
        true
    );
    file_sink->set_level(spdlog::level::debug);

    // Create logger by combining console sink and file sink
    std::vector<spdlog::sink_ptr> sinks{stdout_sink, file_sink};
    auto logger =
        std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);

    // Register logger
    spdlog::set_default_logger(logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");
    spdlog::flush_on(spdlog::level::err); // Ensure error level log is flushed immediately

    spdlog::set_error_handler([](const std::string& errorMessage) {
        spdlog::default_logger()->error(errorMessage);

        // Flush log buffer to ensure all log messages
        // are written to the file before exiting
        spdlog::default_logger()->flush();

        std::exit(EXIT_FAILURE);
    });

    fs::path path = resultDirectory / "time_step_report.csv";
    timeStepReportFile.open(path);
    if (!timeStepReportFile.is_open()) {
        spdlog::error("Could not open the log file: {}", path.string());
    }

    auto writer = csv::make_csv_writer(timeStepReportFile);
    writer << std::vector<std::string>{
        "Time Step",
        "Dt (s)",
        "Time Before (s)",
        "Time After (s)",
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
    if (mps.debugLogCount > 0) {
        spdlog::info(
            "There were {} debug-level logs. Check 'execution.log'.",
            mps.debugLogCount
        );
    } else {
        spdlog::info("There were no debug-level logs.");
    }
    spdlog::info(std::format("Total Simulation Time = {:%T}", totalSimulationTime));

    timeStepReportFile.close();
}

void Simulation::timeStepReport(
    const system_clock::time_point& timeStepStartTime,
    const system_clock::time_point& timeStepEndTime,
    const double& courantNumber,
    const bool& isDebugLogAdded
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

    auto formattedTimeBefore    = std::format("{:.3f}", time - mps.settings.dt);
    auto formattedTimeAfter     = std::format("{:.3f}", time);
    auto formattedElapsed       = std::format("{:%T}", elapsed);
    auto formattedAverage       = std::format("{:.3f}", average);
    auto formattedRemain        = std::format("{:%T}", remain);
    auto formattedLast          = std::format("{:%S}", last);
    auto formattedCourantNumber = std::format("{:.2f}", courantNumber);

    std::cout
        << std::format(
               "{}: dt={}s  t={}s->{}s  fin={}s  elapsed={}  remain={}  ave={}s/step  "
               "last={}s/step  out={}files  Courant={}",
               timeStep,
               mps.settings.dt,
               formattedTimeBefore,
               formattedTimeAfter,
               mps.settings.finishTime,
               formattedElapsed,
               formattedRemain,
               formattedAverage,
               formattedLast,
               outFileNum,
               formattedCourantNumber
           )
        << std::endl;

    auto writer = csv::make_csv_writer(timeStepReportFile);
    writer << std::make_tuple(
        timeStep,
        mps.settings.dt,
        formattedTimeBefore,
        formattedTimeAfter,
        mps.settings.finishTime,
        formattedElapsed,
        elapsed.count(),
        formattedRemain,
        formattedAverage,
        formattedLast,
        outFileNum,
        formattedCourantNumber
    );

    if (isDebugLogAdded) {
        spdlog::debug(
            "Time Step:{} Time:{}s->{}s",
            timeStep,
            formattedTimeBefore,
            formattedTimeAfter
        );
    }
}

bool Simulation::isTimeToExport() {
    double outputInterval = mps.settings.outputInterval;
    double remainder      = std::fmod(time, outputInterval);
    remainder             = std::min(remainder, outputInterval - remainder);
    return remainder < mps.settings.dt / 10.0;
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
