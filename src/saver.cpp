#include "saver.hpp"

#include "csv.hpp"

#include <format>
#include <fstream>

Saver::Saver(const std::filesystem::path& outputDir) {
    this->outputDir = outputDir;
    std::filesystem::create_directory(outputDir / "csv");
}

void Saver::save(
    const std::vector<Particle>& particles, const double& time, const int& fileNum
) {
    std::string filename =
        std::format("{}/csv/output_{:04}.csv", outputDir.string(), fileNum);

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << std::format("Could not open result file: {}", filename) << std::endl;
        exit(-1);
    }
    auto writer = csv::make_csv_writer(outFile);

    writer << std::make_tuple("Time (s)", time);
    writer << std::make_tuple("Number of Particles", particles.size());
    writer << std::make_tuple(
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

        writer << std::make_tuple(
            std::format("{:4}", pi.id),
            std::format("{:2d}", static_cast<int>(pi.type)),
            std::format("{:8.3f}", pi.position.x()),
            std::format("{:8.3f}", pi.position.y()),
            std::format("{:8.3f}", pi.position.z()),
            std::format("{:8.3f}", pi.velocity.x()),
            std::format("{:8.3f}", pi.velocity.y()),
            std::format("{:8.3f}", pi.velocity.z()),
            std::format("{:9.3f}", pi.pressure),
            std::format("{:8.3f}", pi.numberDensity)
        );
    }

    outFile.close();
}
