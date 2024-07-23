#include "saver.hpp"

#include "csv.hpp"

#include <format>
#include <fstream>

using std::endl;
using std::format;
using std::make_tuple;

Saver::Saver(const std::filesystem::path& outputDir) {
    this->outputDir = outputDir;
    std::filesystem::create_directory(outputDir / "csv");
    std::filesystem::create_directory(outputDir / "vtu");
}

void Saver::save(const std::vector<Particle>& particles, const double& time) {
    toCsv(particles, time);
    toVtu(particles, time);
    numberOfFiles++;
}

void Saver::toCsv(const std::vector<Particle>& particles, const double& time) {
    std::string filename =
        format("{}/csv/output_{:04}.csv", outputDir.string(), numberOfFiles);

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
        "Number Density Ratio"
    );
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost) {
            continue;
        }

        writer << make_tuple(
            pi.id,
            static_cast<int>(pi.type),
            pi.position.x(),
            pi.position.y(),
            pi.position.z(),
            pi.velocity.x(),
            pi.velocity.y(),
            pi.velocity.z(),
            pi.pressure,
            pi.numberDensityRatio
        );
    }

    outFile.close();
}

void Saver::toVtu(const std::vector<Particle>& particles, const double& time) {
    std::string filename =
        format("{}/vtu/output_{:04}.vtu", outputDir.string(), numberOfFiles);

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << format("Could not open result file: {}", filename) << std::endl;
        exit(-1);
    }

    // --------------
    // --- Header ---
    // --------------
    outFile << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
    outFile
        << "<VTKFile byte_order='LittleEndian' version='0.1' type = 'UnstructuredGrid'>"
        << endl;
    outFile << "<UnstructuredGrid>" << endl;

    // -----------------
    // --- Time data ---
    // -----------------
    outFile << "<FieldData>" << endl;
    outFile << "<DataArray type='Float64' Name='Time' NumberOfComponents='1' "
               "NumberOfTuples='1' format='ascii'>"
            << endl;
    outFile << time << endl;
    outFile << "</DataArray>" << endl;
    outFile << "</FieldData>" << endl;

    /// ------------------
    /// ----- Points -----
    /// ------------------
    outFile << "<Piece NumberOfCells='" << particles.size() << "' NumberOfPoints='"
            << particles.size() << "'>" << endl;
    outFile << "<Points>" << endl;
    outFile << "<DataArray NumberOfComponents='3' type='Float64' "
               "Name='position' format='ascii'>"
            << endl;
    for (const auto& pi : particles) {
        outFile << pi.position.x() << " ";
        outFile << pi.position.y() << " ";
        outFile << pi.position.z() << endl;
    }
    outFile << "</DataArray>" << endl;
    outFile << "</Points>" << endl;

    // ---------------------
    // ----- PointData -----
    // ---------------------
    outFile << "<PointData>" << endl;

    dataArrayBegin(outFile, "1", "Int32", "Type");
    for (auto& pi : particles) {
        outFile << static_cast<int>(pi.type) << endl;
    }
    dataArrayEnd(outFile);

    dataArrayBegin(outFile, "3", "Float64", "Velocity");
    for (const auto& pi : particles) {
        outFile << pi.velocity.x() << " ";
        outFile << pi.velocity.y() << " ";
        outFile << 0.0 << endl;
    }
    dataArrayEnd(outFile);

    dataArrayBegin(outFile, "1", "Float64", "Pressure");
    for (const auto& pi : particles)
        outFile << pi.pressure << endl;
    dataArrayEnd(outFile);

    dataArrayBegin(outFile, "1", "Float64", "Number Density Ratio");
    for (auto& pi : particles)
        outFile << pi.numberDensityRatio << std::endl;
    dataArrayEnd(outFile);

    dataArrayBegin(outFile, "1", "Int32", "Boundary Condition");
    for (auto& pi : particles)
        outFile << static_cast<int>(pi.boundaryCondition) << std::endl;
    dataArrayEnd(outFile);

    outFile << "</PointData>" << endl;

    // -----------------
    // ----- Cells -----
    // -----------------
    outFile << "<Cells>" << endl;
    outFile << "<DataArray type='Int32' Name='connectivity' format='ascii'>" << endl;
    for (int i = 0; i < particles.size(); i++) {
        outFile << i << endl;
    }
    outFile << "</DataArray>" << endl;
    outFile << "<DataArray type='Int32' Name='offsets' format='ascii'>" << endl;
    for (int i = 0; i < particles.size(); i++) {
        outFile << i + 1 << endl;
    }
    outFile << "</DataArray>" << endl;
    outFile << "<DataArray type='UInt8' Name='types' format='ascii'>" << endl;
    for (int i = 0; i < particles.size(); i++) {
        outFile << "1" << endl;
    }
    outFile << "</DataArray>" << endl;
    outFile << "</Cells>" << endl;

    // ------------------
    // ----- Footer -----
    // ------------------
    outFile << "</Piece>" << endl;
    outFile << "</UnstructuredGrid>" << endl;
    outFile << "</VTKFile>" << endl;
}

void Saver::dataArrayBegin(
    std::ofstream& outFile,
    const std::string& numberOfComponents,
    const std::string& type,
    const std::string& name
) {
    outFile << "<DataArray NumberOfComponents='" << numberOfComponents << "' type='"
            << type << "' Name='" << name << "' format='ascii'>" << endl;
}

void Saver::dataArrayEnd(std::ofstream& outFile) {
    outFile << "</DataArray>" << endl;
}
