#pragma once
#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include <cstdint>
#include <string>
#include <vector>

struct CompressionStats {
    int original_record_length = 0;
    int compressed_record_length = 0;
    double original_cutoff_energy = 0.0;
    double compressed_cutoff_energy = 0.0;
    std::uint64_t original_file_size = 0;
    std::uint64_t compressed_file_size = 0;
    std::uint64_t original_coefficients = 0;
    std::uint64_t kept_coefficients = 0;
    std::uint64_t zeroed_components = 0;
};

int compressWavecar(
    const std::string& input_path,
    const std::string& output_path,
    double cutoff_fraction,
    bool zero_small_values,
    CompressionStats* stats = nullptr
);

struct AnalysisPoint {
    double cutoff_fraction = 0.0;
    double cutoff_energy = 0.0;
    std::uint64_t kept_coefficients = 0;
    std::uint64_t total_coefficients = 0;
    double coefficient_fraction = 0.0;
    double retained_norm_fraction = 0.0;
    double lost_norm_fraction = 0.0;
    double min_band_retained_norm_fraction = 0.0;
    int estimated_record_length = 0;
    std::uint64_t estimated_file_size = 0;
    double estimated_size_ratio = 0.0;
};

int analyzeWavecar(
    const std::string& input_path,
    const std::vector<double>& cutoff_fractions,
    std::vector<AnalysisPoint>& analysis
);

#endif // COMPRESSOR_HPP
