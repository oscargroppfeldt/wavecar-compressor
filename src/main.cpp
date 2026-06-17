#include "compressor.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void printUsage(const char* program_name) {
    std::cerr << "Usage:\n"
              << "  " << program_name
              << " compress <input WAVECAR> <output WAVECAR> <cutoff_fraction> [zero_small]\n\n"
              << "  " << program_name
              << " analyze <input WAVECAR> [cutoff_fraction ...]\n\n"
              << "Arguments:\n"
              << "  cutoff_fraction  value in [0, 1]; keeps plane waves with E <= value * ENCUT\n"
              << "  zero_small       optional boolean: true/false, 1/0, yes/no; default false\n"
              << "  analyze          prints CSV with retained coefficient norm and size estimates\n";
}

bool parseBool(const std::string& value, bool& parsed_value) {
    if (value == "1" || value == "true" || value == "TRUE" || value == "yes" ||
        value == "YES" || value == "on" || value == "ON") {
        parsed_value = true;
        return true;
    }

    if (value == "0" || value == "false" || value == "FALSE" || value == "no" ||
        value == "NO" || value == "off" || value == "OFF") {
        parsed_value = false;
        return true;
    }

    return false;
}

bool parseDouble(const std::string& value, double& parsed_value) {
    try {
        std::size_t parsed_chars = 0;
        parsed_value = std::stod(value, &parsed_chars);
        return parsed_chars == value.size();
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    }
}

std::vector<double> defaultAnalysisFractions() {
    return {0.0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99, 1.0};
}

int runCompress(int argc, char* argv[]) {
    if (argc != 5 && argc != 6) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const std::string input_path = argv[2];
    const std::string output_path = argv[3];

    double cutoff_fraction = 0.0;
    if (!parseDouble(argv[4], cutoff_fraction)) {
        std::cerr << "Invalid cutoff fraction: " << argv[4] << std::endl;
        return EXIT_FAILURE;
    }

    bool zero_small_values = false;
    if (argc == 6 && !parseBool(argv[5], zero_small_values)) {
        std::cerr << "Invalid zero_small value: " << argv[5] << std::endl;
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    CompressionStats stats;
    const int status = compressWavecar(
        input_path,
        output_path,
        cutoff_fraction,
        zero_small_values,
        &stats
    );

    if (status != 0) {
        return EXIT_FAILURE;
    }

    std::cout << "Compression complete\n"
              << "  ENCUT: " << stats.original_cutoff_energy << " -> "
              << stats.compressed_cutoff_energy << "\n"
              << "  record length: " << stats.original_record_length << " -> "
              << stats.compressed_record_length << "\n"
              << "  coefficients: " << stats.kept_coefficients << " / "
              << stats.original_coefficients << " kept\n"
              << "  file size: " << stats.original_file_size << " -> "
              << stats.compressed_file_size << " bytes\n"
              << "  compression ratio: " << static_cast<double>(stats.compressed_file_size) /
              static_cast<double>(stats.original_file_size) << "\n";

    if (zero_small_values) {
        std::cout << "  zeroed components: " << stats.zeroed_components << "\n";
    }

    return EXIT_SUCCESS;
}

int runAnalyze(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    std::vector<double> cutoff_fractions;
    if (argc == 3) {
        cutoff_fractions = defaultAnalysisFractions();
    } else {
        for (int idx = 3; idx < argc; ++idx) {
            double cutoff_fraction = 0.0;
            if (!parseDouble(argv[idx], cutoff_fraction) || !std::isfinite(cutoff_fraction)) {
                std::cerr << "Invalid cutoff fraction: " << argv[idx] << std::endl;
                return EXIT_FAILURE;
            }
            cutoff_fractions.push_back(cutoff_fraction);
        }
    }

    std::vector<AnalysisPoint> analysis;
    const int status = analyzeWavecar(argv[2], cutoff_fractions, analysis);
    if (status != 0) {
        return EXIT_FAILURE;
    }

    std::cout << std::setprecision(17)
              << "cutoff_fraction,encut,kept_coefficients,total_coefficients,"
              << "coefficient_fraction,retained_norm_fraction,lost_norm_fraction,"
              << "min_band_retained_norm_fraction,estimated_record_length,"
              << "estimated_file_size,estimated_size_ratio\n";

    for (const AnalysisPoint& point : analysis) {
        std::cout << point.cutoff_fraction << ','
                  << point.cutoff_energy << ','
                  << point.kept_coefficients << ','
                  << point.total_coefficients << ','
                  << point.coefficient_fraction << ','
                  << point.retained_norm_fraction << ','
                  << point.lost_norm_fraction << ','
                  << point.min_band_retained_norm_fraction << ','
                  << point.estimated_record_length << ','
                  << point.estimated_file_size << ','
                  << point.estimated_size_ratio << '\n';
    }

    return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const std::string command = argv[1];
    if (command == "compress") {
        return runCompress(argc, argv);
    }
    if (command == "analyze") {
        return runAnalyze(argc, argv);
    }

    printUsage(argv[0]);
    return EXIT_FAILURE;
}
