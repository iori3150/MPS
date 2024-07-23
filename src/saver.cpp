#include "saver.hpp"

#include "csv.hpp"

#include <format>
#include <fstream>

using std::format;
using std::make_tuple;

Saver::Saver(const std::filesystem::path& outputDir) {
    this->outputDir = outputDir;
    std::filesystem::create_directory(outputDir / "csv");
}

void Saver::save(
    const std::vector<Particle>& particles, const double& time, const int& fileNum
) {
    std::string filename = format("{}/csv/output_{:04}.csv", outputDir.string(), fileNum);

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << format("Could not open result file: {}", filename) << std::endl;
        exit(-1);
    }
    auto writer = csv::make_csv_writer(outFile);

    writer << make_tuple("Time (s)", time);
    writer << make_tuple("Number of Particles", particles.size());
    writer << make_tuple(
        "ID",
        "Type",
        "Position.x (m)",
        "Position.y (m)",
        "Position.z (m)",
        "Velocity.x (m/s)",
        "Velocity.y (m/s)",
        "Velocity.z (m/s)",
        "Pressure (Pa)",
        "Number Density"
    );
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost) {
            continue;
        }

        writer << make_tuple(
            format("{:4}", pi.id),
            format("{:2d}", static_cast<int>(pi.type)),
            format("{:8.3f}", pi.position.x()),
            format("{:8.3f}", pi.position.y()),
            format("{:8.3f}", pi.position.z()),
            format("{:8.3f}", pi.velocity.x()),
            format("{:8.3f}", pi.velocity.y()),
            format("{:8.3f}", pi.velocity.z()),
            format("{:9.3f}", pi.pressure),
            format("{:8.3f}", pi.numberDensity)
        );
    }

    outFile.close();
}
