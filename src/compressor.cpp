#include "compressor.hpp"
#include "reader.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace {

constexpr double kHbarSquaredOver2M = 3.80998212; // eV Angstrom^2
constexpr double kTwoPi = 6.28318530717958647692;
constexpr double kKpointTolerance = 1.0e-12;

struct CoefficientLayout {
    std::size_t component_bytes = 0;
    std::size_t complex_bytes = 0;
};

enum class StorageMode {
    Full,
    GammaHalfZ,
    GammaHalfX,
    GammaHalfY,
};

struct PlaneWaveMask {
    Vec3 kpoint;
    int original_count = 0;
    int kept_count = 0;
    StorageMode storage_mode = StorageMode::Full;
    std::vector<unsigned char> keep;
};

struct AnalysisPlaneWaveMap {
    Vec3 kpoint;
    int original_count = 0;
    StorageMode storage_mode = StorageMode::Full;
    std::vector<std::size_t> first_retained_fraction;
    std::vector<int> kept_counts;
};

bool coefficientLayout(int rtag, CoefficientLayout& layout) {
    switch (rtag) {
        case 45200:
            layout = {sizeof(float), 2 * sizeof(float)};
            return true;
        case 45210:
            layout = {sizeof(double), 2 * sizeof(double)};
            return true;
        default:
            return false;
    }
}

double dot(const Vec3& left, const Vec3& right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

Vec3 cross(const Vec3& left, const Vec3& right) {
    return Vec3(
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x
    );
}

Vec3 scale(const Vec3& vector, double scalar) {
    return Vec3(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

double norm(const Vec3& vector) {
    return std::sqrt(dot(vector, vector));
}

std::array<Vec3, 3> reciprocalLattice(const Cell& lattice) {
    const double volume = dot(lattice.a, cross(lattice.b, lattice.c));
    const double factor = kTwoPi / volume;
    return {
        scale(cross(lattice.b, lattice.c), factor),
        scale(cross(lattice.c, lattice.a), factor),
        scale(cross(lattice.a, lattice.b), factor),
    };
}

bool isGammaKpoint(const Vec3& kpoint) {
    return std::abs(kpoint.x) < kKpointTolerance &&
           std::abs(kpoint.y) < kKpointTolerance &&
           std::abs(kpoint.z) < kKpointTolerance;
}

bool withinCutoff(double energy, double cutoff) {
    const double tolerance = std::max(1.0, cutoff) * 1.0e-10;
    return energy <= cutoff + tolerance;
}

int wrapIndex(int raw_index, int max_index) {
    if (raw_index <= max_index) {
        return raw_index;
    }
    return raw_index - 2 * max_index - 1;
}

bool storesPlaneWave(StorageMode mode, int g1, int g2, int g3) {
    switch (mode) {
        case StorageMode::Full:
            return true;
        case StorageMode::GammaHalfZ:
            return g3 > 0 || (g3 == 0 && (g2 > 0 || (g2 == 0 && g1 >= 0)));
        case StorageMode::GammaHalfX:
            return g1 > 0 || (g1 == 0 && (g2 > 0 || (g2 == 0 && g3 >= 0)));
        case StorageMode::GammaHalfY:
            return g2 > 0 || (g2 == 0 && (g3 > 0 || (g3 == 0 && g1 >= 0)));
    }
    return false;
}

const char* storageModeName(StorageMode mode) {
    switch (mode) {
        case StorageMode::Full:
            return "full";
        case StorageMode::GammaHalfZ:
            return "gamma-half-z";
        case StorageMode::GammaHalfX:
            return "gamma-half-x";
        case StorageMode::GammaHalfY:
            return "gamma-half-y";
    }
    return "unknown";
}

std::vector<double> planeWaveEnergies(
    const WavecarMetadata& metadata,
    const Vec3& kpoint,
    double cutoff,
    StorageMode storage_mode
) {
    const std::array<Vec3, 3> reciprocal = reciprocalLattice(metadata.lattice);
    const double wave_number_cutoff = std::sqrt(cutoff / kHbarSquaredOver2M);
    const std::array<int, 3> max_indices = {
        static_cast<int>(wave_number_cutoff / norm(reciprocal[0])) + 1,
        static_cast<int>(wave_number_cutoff / norm(reciprocal[1])) + 1,
        static_cast<int>(wave_number_cutoff / norm(reciprocal[2])) + 1,
    };

    std::vector<double> energies;

    for (int raw_g3 = 0; raw_g3 <= 2 * max_indices[2]; ++raw_g3) {
        const int g3 = wrapIndex(raw_g3, max_indices[2]);
        for (int raw_g2 = 0; raw_g2 <= 2 * max_indices[1]; ++raw_g2) {
            const int g2 = wrapIndex(raw_g2, max_indices[1]);
            for (int raw_g1 = 0; raw_g1 <= 2 * max_indices[0]; ++raw_g1) {
                const int g1 = wrapIndex(raw_g1, max_indices[0]);
                if (!storesPlaneWave(storage_mode, g1, g2, g3)) {
                    continue;
                }

                const Vec3 wave_vector(
                    (g1 + kpoint.x) * reciprocal[0].x +
                        (g2 + kpoint.y) * reciprocal[1].x +
                        (g3 + kpoint.z) * reciprocal[2].x,
                    (g1 + kpoint.x) * reciprocal[0].y +
                        (g2 + kpoint.y) * reciprocal[1].y +
                        (g3 + kpoint.z) * reciprocal[2].y,
                    (g1 + kpoint.x) * reciprocal[0].z +
                        (g2 + kpoint.y) * reciprocal[1].z +
                        (g3 + kpoint.z) * reciprocal[2].z
                );
                const double energy = kHbarSquaredOver2M * dot(wave_vector, wave_vector);
                if (withinCutoff(energy, cutoff)) {
                    energies.push_back(energy);
                }
            }
        }
    }

    return energies;
}

bool matchingPlaneWaves(
    const WavecarMetadata& metadata,
    const Vec3& kpoint,
    int stored_count,
    StorageMode& storage_mode,
    std::vector<double>& energies
) {
    energies = planeWaveEnergies(metadata, kpoint, metadata.cutoff_energy, StorageMode::Full);
    if (static_cast<int>(energies.size()) == stored_count) {
        storage_mode = StorageMode::Full;
        return true;
    }

    if (!isGammaKpoint(kpoint)) {
        std::cerr << "Plane-wave count mismatch: file has " << stored_count
                  << ", full-grid enumeration found " << energies.size() << std::endl;
        return false;
    }

    const std::array<StorageMode, 3> gamma_modes = {
        StorageMode::GammaHalfZ,
        StorageMode::GammaHalfX,
        StorageMode::GammaHalfY,
    };

    for (StorageMode mode : gamma_modes) {
        energies = planeWaveEnergies(metadata, kpoint, metadata.cutoff_energy, mode);
        if (static_cast<int>(energies.size()) == stored_count) {
            storage_mode = mode;
            return true;
        }
    }

    std::cerr << "Plane-wave count mismatch: file has " << stored_count
              << ", and neither full nor gamma-half enumeration matched." << std::endl;
    return false;
}

std::uint64_t kpointRecordIndex(
    const WavecarMetadata& metadata,
    int spin_idx,
    int kpoint_idx
) {
    return 2ULL +
           static_cast<std::uint64_t>(spin_idx * metadata.kpoint_count + kpoint_idx) *
               static_cast<std::uint64_t>(metadata.band_count + 1);
}

std::streamoff recordOffset(std::uint64_t record_index, int record_length) {
    return static_cast<std::streamoff>(record_index) *
           static_cast<std::streamoff>(record_length);
}

bool readAt(std::ifstream& input, std::streamoff offset, char* data, std::size_t bytes) {
    input.seekg(offset, std::ios::beg);
    if (!input) {
        return false;
    }
    input.read(data, static_cast<std::streamsize>(bytes));
    return static_cast<bool>(input);
}

void putDouble(std::vector<char>& record, std::size_t double_index, double value) {
    std::memcpy(record.data() + double_index * sizeof(double), &value, sizeof(double));
}

void copyCoefficient(
    const char* source,
    char* destination,
    const CoefficientLayout& layout,
    bool zero_small_values,
    std::uint64_t& zeroed_components
) {
    if (!zero_small_values) {
        std::memcpy(destination, source, layout.complex_bytes);
        return;
    }

    const double threshold = std::numeric_limits<double>::epsilon();
    if (layout.component_bytes == sizeof(float)) {
        float values[2] = {0.0F, 0.0F};
        std::memcpy(values, source, sizeof(values));
        for (float& value : values) {
            if (value != 0.0F && std::abs(static_cast<double>(value)) < threshold) {
                value = 0.0F;
                ++zeroed_components;
            }
        }
        std::memcpy(destination, values, sizeof(values));
        return;
    }

    double values[2] = {0.0, 0.0};
    std::memcpy(values, source, sizeof(values));
    for (double& value : values) {
        if (value != 0.0 && std::abs(value) < threshold) {
            value = 0.0;
            ++zeroed_components;
        }
    }
    std::memcpy(destination, values, sizeof(values));
}

double coefficientNormSquared(const char* source, const CoefficientLayout& layout) {
    if (layout.component_bytes == sizeof(float)) {
        float values[2] = {0.0F, 0.0F};
        std::memcpy(values, source, sizeof(values));
        const double real = static_cast<double>(values[0]);
        const double imag = static_cast<double>(values[1]);
        return real * real + imag * imag;
    }

    double values[2] = {0.0, 0.0};
    std::memcpy(values, source, sizeof(values));
    return values[0] * values[0] + values[1] * values[1];
}

bool parseKpointHeader(
    std::ifstream& input,
    const WavecarMetadata& metadata,
    int spin_idx,
    int kpoint_idx,
    int& coefficient_count,
    Vec3& kpoint
) {
    double values[4] = {0.0, 0.0, 0.0, 0.0};
    const std::streamoff offset = recordOffset(
        kpointRecordIndex(metadata, spin_idx, kpoint_idx),
        metadata.record_length
    );
    if (!readAt(input, offset, reinterpret_cast<char*>(values), sizeof(values))) {
        return false;
    }

    coefficient_count = static_cast<int>(values[0]);
    kpoint = Vec3(values[1], values[2], values[3]);
    return coefficient_count > 0;
}

bool validateCutoffFractions(const std::vector<double>& cutoff_fractions) {
    if (cutoff_fractions.empty()) {
        std::cerr << "At least one cutoff fraction is required." << std::endl;
        return false;
    }

    for (double fraction : cutoff_fractions) {
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
            std::cerr << "Cutoff fractions must be finite values in [0, 1]." << std::endl;
            return false;
        }
    }

    return true;
}

std::vector<double> sortedUniqueFractions(const std::vector<double>& cutoff_fractions) {
    std::vector<double> fractions = cutoff_fractions;
    std::sort(fractions.begin(), fractions.end());
    fractions.erase(
        std::unique(
            fractions.begin(),
            fractions.end(),
            [](double left, double right) {
                return std::abs(left - right) < 1.0e-14;
            }
        ),
        fractions.end()
    );
    return fractions;
}

} // namespace

int compressWavecar(
    const std::string& input_path,
    const std::string& output_path,
    double cutoff_fraction,
    bool zero_small_values,
    CompressionStats* stats
) {
    if (!std::isfinite(cutoff_fraction) || cutoff_fraction < 0.0 || cutoff_fraction > 1.0) {
        std::cerr << "Cutoff fraction must be in the range [0, 1]." << std::endl;
        return -1;
    }

    if (input_path == output_path) {
        std::cerr << "Input and output paths must be different." << std::endl;
        return -2;
    }

    namespace fs = std::filesystem;
    std::error_code fs_error;
    if (fs::exists(output_path, fs_error) &&
        fs::equivalent(input_path, output_path, fs_error)) {
        std::cerr << "Input and output paths refer to the same file." << std::endl;
        return -2;
    }

    std::ifstream input(input_path, std::ios::binary | std::ios::in);
    if (!input) {
        std::cerr << "Error opening input file: " << input_path << std::endl;
        return -3;
    }

    WavecarMetadata metadata;
    const int metadata_status = readMetadata(metadata, input);
    if (metadata_status != 0) {
        std::cerr << "Error reading WAVECAR metadata from: " << input_path << std::endl;
        return metadata_status;
    }

    CoefficientLayout coefficient_layout;
    if (!coefficientLayout(metadata.RTAG_TAG, coefficient_layout)) {
        std::cerr << "Unsupported WAVECAR RTAG: " << metadata.RTAG_TAG << std::endl;
        return -4;
    }

    const double target_cutoff = cutoff_fraction * metadata.cutoff_energy;
    const std::size_t kpoint_header_bytes =
        static_cast<std::size_t>(4 + 3 * metadata.band_count) * sizeof(double);
    const std::size_t mask_count =
        static_cast<std::size_t>(metadata.spin_count * metadata.kpoint_count);
    std::vector<PlaneWaveMask> masks(mask_count);

    int max_kept_coefficients = 0;
    for (int spin_idx = 0; spin_idx < metadata.spin_count; ++spin_idx) {
        for (int kpoint_idx = 0; kpoint_idx < metadata.kpoint_count; ++kpoint_idx) {
            PlaneWaveMask& mask =
                masks[static_cast<std::size_t>(spin_idx * metadata.kpoint_count + kpoint_idx)];

            if (!parseKpointHeader(
                    input,
                    metadata,
                    spin_idx,
                    kpoint_idx,
                    mask.original_count,
                    mask.kpoint
                )) {
                std::cerr << "Error reading k-point header for spin " << spin_idx
                          << ", k-point " << kpoint_idx << std::endl;
                return -5;
            }

            std::vector<double> energies;
            if (!matchingPlaneWaves(
                    metadata,
                    mask.kpoint,
                    mask.original_count,
                    mask.storage_mode,
                    energies
                )) {
                std::cerr << "Unable to build plane-wave mask for spin " << spin_idx
                          << ", k-point " << kpoint_idx << std::endl;
                return -6;
            }

            mask.keep.resize(energies.size(), 0U);
            for (std::size_t idx = 0; idx < energies.size(); ++idx) {
                if (withinCutoff(energies[idx], target_cutoff)) {
                    mask.keep[idx] = 1U;
                    ++mask.kept_count;
                }
            }
            max_kept_coefficients = std::max(max_kept_coefficients, mask.kept_count);
        }
    }

    const std::size_t coefficient_record_bytes =
        static_cast<std::size_t>(max_kept_coefficients + 1) *
        coefficient_layout.complex_bytes;
    const int compressed_record_length = static_cast<int>(
        std::max(kpoint_header_bytes, coefficient_record_bytes)
    );

    std::ofstream output(output_path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!output) {
        std::cerr << "Error creating output file: " << output_path << std::endl;
        return -7;
    }

    std::vector<char> record(static_cast<std::size_t>(compressed_record_length), 0);
    putDouble(record, 0, static_cast<double>(compressed_record_length));
    putDouble(record, 1, static_cast<double>(metadata.spin_count));
    putDouble(record, 2, static_cast<double>(metadata.RTAG_TAG));
    output.write(record.data(), static_cast<std::streamsize>(record.size()));

    std::fill(record.begin(), record.end(), 0);
    putDouble(record, 0, static_cast<double>(metadata.kpoint_count));
    putDouble(record, 1, static_cast<double>(metadata.band_count));
    putDouble(record, 2, target_cutoff);
    putDouble(record, 3, metadata.lattice.a.x);
    putDouble(record, 4, metadata.lattice.a.y);
    putDouble(record, 5, metadata.lattice.a.z);
    putDouble(record, 6, metadata.lattice.b.x);
    putDouble(record, 7, metadata.lattice.b.y);
    putDouble(record, 8, metadata.lattice.b.z);
    putDouble(record, 9, metadata.lattice.c.x);
    putDouble(record, 10, metadata.lattice.c.y);
    putDouble(record, 11, metadata.lattice.c.z);
    putDouble(record, 12, metadata.fermi_energy);
    output.write(record.data(), static_cast<std::streamsize>(record.size()));

    if (!output) {
        std::cerr << "Error writing WAVECAR header to: " << output_path << std::endl;
        return -8;
    }

    CompressionStats local_stats;
    local_stats.original_record_length = metadata.record_length;
    local_stats.compressed_record_length = compressed_record_length;
    local_stats.original_cutoff_energy = metadata.cutoff_energy;
    local_stats.compressed_cutoff_energy = target_cutoff;
    local_stats.original_file_size = fs::file_size(input_path, fs_error);

    for (int spin_idx = 0; spin_idx < metadata.spin_count; ++spin_idx) {
        for (int kpoint_idx = 0; kpoint_idx < metadata.kpoint_count; ++kpoint_idx) {
            const std::size_t mask_idx =
                static_cast<std::size_t>(spin_idx * metadata.kpoint_count + kpoint_idx);
            const PlaneWaveMask& mask = masks[mask_idx];
            const std::uint64_t input_kpoint_record =
                kpointRecordIndex(metadata, spin_idx, kpoint_idx);

            std::fill(record.begin(), record.end(), 0);
            if (!readAt(
                    input,
                    recordOffset(input_kpoint_record, metadata.record_length),
                    record.data(),
                    kpoint_header_bytes
                )) {
                std::cerr << "Error reading k-point payload for spin " << spin_idx
                          << ", k-point " << kpoint_idx << std::endl;
                return -9;
            }
            putDouble(record, 0, static_cast<double>(mask.kept_count));
            output.write(record.data(), static_cast<std::streamsize>(record.size()));

            const std::size_t input_coefficient_bytes =
                static_cast<std::size_t>(mask.original_count) *
                coefficient_layout.complex_bytes;
            if (input_coefficient_bytes >
                static_cast<std::size_t>(metadata.record_length)) {
                std::cerr << "Coefficient record is larger than WAVECAR record length."
                          << std::endl;
                return -10;
            }

            std::vector<char> input_coefficients(input_coefficient_bytes, 0);
            std::vector<char> output_coefficients(
                static_cast<std::size_t>(compressed_record_length),
                0
            );

            for (int band_idx = 0; band_idx < metadata.band_count; ++band_idx) {
                const std::uint64_t input_band_record =
                    input_kpoint_record + 1ULL + static_cast<std::uint64_t>(band_idx);
                if (!readAt(
                        input,
                        recordOffset(input_band_record, metadata.record_length),
                        input_coefficients.data(),
                        input_coefficients.size()
                    )) {
                    std::cerr << "Error reading coefficients for spin " << spin_idx
                              << ", k-point " << kpoint_idx << ", band " << band_idx
                              << std::endl;
                    return -11;
                }

                std::fill(output_coefficients.begin(), output_coefficients.end(), 0);
                int output_coefficient_idx = 0;
                for (int coeff_idx = 0; coeff_idx < mask.original_count; ++coeff_idx) {
                    if (mask.keep[static_cast<std::size_t>(coeff_idx)] == 0U) {
                        continue;
                    }

                    const std::size_t source_offset =
                        static_cast<std::size_t>(coeff_idx) *
                        coefficient_layout.complex_bytes;
                    const std::size_t destination_offset =
                        static_cast<std::size_t>(output_coefficient_idx) *
                        coefficient_layout.complex_bytes;
                    copyCoefficient(
                        input_coefficients.data() + source_offset,
                        output_coefficients.data() + destination_offset,
                        coefficient_layout,
                        zero_small_values,
                        local_stats.zeroed_components
                    );
                    ++output_coefficient_idx;
                }

                output.write(
                    output_coefficients.data(),
                    static_cast<std::streamsize>(output_coefficients.size())
                );
                if (!output) {
                    std::cerr << "Error writing coefficients for spin " << spin_idx
                              << ", k-point " << kpoint_idx << ", band " << band_idx
                              << std::endl;
                    return -12;
                }

                local_stats.original_coefficients +=
                    static_cast<std::uint64_t>(mask.original_count);
                local_stats.kept_coefficients +=
                    static_cast<std::uint64_t>(mask.kept_count);
            }

            std::cout << "spin " << spin_idx << ", k-point " << kpoint_idx
                      << ": kept " << mask.kept_count << " / " << mask.original_count
                      << " plane waves (" << storageModeName(mask.storage_mode) << ")"
                      << std::endl;
        }
    }

    output.close();
    if (!output) {
        std::cerr << "Error finalizing output file: " << output_path << std::endl;
        return -13;
    }

    local_stats.compressed_file_size = fs::file_size(output_path, fs_error);
    if (stats != nullptr) {
        *stats = local_stats;
    }

    return 0;
}

int analyzeWavecar(
    const std::string& input_path,
    const std::vector<double>& cutoff_fractions,
    std::vector<AnalysisPoint>& analysis
) {
    analysis.clear();
    if (!validateCutoffFractions(cutoff_fractions)) {
        return -1;
    }

    const std::vector<double> fractions = sortedUniqueFractions(cutoff_fractions);

    std::ifstream input(input_path, std::ios::binary | std::ios::in);
    if (!input) {
        std::cerr << "Error opening input file: " << input_path << std::endl;
        return -2;
    }

    WavecarMetadata metadata;
    const int metadata_status = readMetadata(metadata, input);
    if (metadata_status != 0) {
        std::cerr << "Error reading WAVECAR metadata from: " << input_path << std::endl;
        return metadata_status;
    }

    CoefficientLayout coefficient_layout;
    if (!coefficientLayout(metadata.RTAG_TAG, coefficient_layout)) {
        std::cerr << "Unsupported WAVECAR RTAG: " << metadata.RTAG_TAG << std::endl;
        return -3;
    }

    const std::size_t fraction_count = fractions.size();
    const std::size_t map_count =
        static_cast<std::size_t>(metadata.spin_count * metadata.kpoint_count);
    std::vector<AnalysisPlaneWaveMap> maps(map_count);
    std::vector<int> max_kept_coefficients(fraction_count, 0);
    std::uint64_t total_coefficients = 0;

    for (int spin_idx = 0; spin_idx < metadata.spin_count; ++spin_idx) {
        for (int kpoint_idx = 0; kpoint_idx < metadata.kpoint_count; ++kpoint_idx) {
            AnalysisPlaneWaveMap& map =
                maps[static_cast<std::size_t>(spin_idx * metadata.kpoint_count + kpoint_idx)];

            if (!parseKpointHeader(
                    input,
                    metadata,
                    spin_idx,
                    kpoint_idx,
                    map.original_count,
                    map.kpoint
                )) {
                std::cerr << "Error reading k-point header for spin " << spin_idx
                          << ", k-point " << kpoint_idx << std::endl;
                return -4;
            }

            std::vector<double> energies;
            if (!matchingPlaneWaves(
                    metadata,
                    map.kpoint,
                    map.original_count,
                    map.storage_mode,
                    energies
                )) {
                std::cerr << "Unable to build plane-wave map for spin " << spin_idx
                          << ", k-point " << kpoint_idx << std::endl;
                return -5;
            }

            map.first_retained_fraction.resize(energies.size(), fraction_count);
            std::vector<int> count_diff(fraction_count + 1, 0);

            for (std::size_t coeff_idx = 0; coeff_idx < energies.size(); ++coeff_idx) {
                for (std::size_t fraction_idx = 0; fraction_idx < fraction_count; ++fraction_idx) {
                    const double target_cutoff =
                        fractions[fraction_idx] * metadata.cutoff_energy;
                    if (withinCutoff(energies[coeff_idx], target_cutoff)) {
                        map.first_retained_fraction[coeff_idx] = fraction_idx;
                        ++count_diff[fraction_idx];
                        break;
                    }
                }
            }

            map.kept_counts.assign(fraction_count, 0);
            int running_count = 0;
            for (std::size_t fraction_idx = 0; fraction_idx < fraction_count; ++fraction_idx) {
                running_count += count_diff[fraction_idx];
                map.kept_counts[fraction_idx] = running_count;
                max_kept_coefficients[fraction_idx] = std::max(
                    max_kept_coefficients[fraction_idx],
                    running_count
                );
            }

            total_coefficients +=
                static_cast<std::uint64_t>(map.original_count) *
                static_cast<std::uint64_t>(metadata.band_count);
        }
    }

    std::vector<long double> retained_norm(fraction_count, 0.0L);
    std::vector<double> min_band_retained_norm(fraction_count, 1.0);
    long double total_norm = 0.0L;

    for (int spin_idx = 0; spin_idx < metadata.spin_count; ++spin_idx) {
        for (int kpoint_idx = 0; kpoint_idx < metadata.kpoint_count; ++kpoint_idx) {
            const std::size_t map_idx =
                static_cast<std::size_t>(spin_idx * metadata.kpoint_count + kpoint_idx);
            const AnalysisPlaneWaveMap& map = maps[map_idx];
            const std::uint64_t input_kpoint_record =
                kpointRecordIndex(metadata, spin_idx, kpoint_idx);

            const std::size_t input_coefficient_bytes =
                static_cast<std::size_t>(map.original_count) *
                coefficient_layout.complex_bytes;
            if (input_coefficient_bytes >
                static_cast<std::size_t>(metadata.record_length)) {
                std::cerr << "Coefficient record is larger than WAVECAR record length."
                          << std::endl;
                return -6;
            }

            std::vector<char> input_coefficients(input_coefficient_bytes, 0);
            for (int band_idx = 0; band_idx < metadata.band_count; ++band_idx) {
                const std::uint64_t input_band_record =
                    input_kpoint_record + 1ULL + static_cast<std::uint64_t>(band_idx);
                if (!readAt(
                        input,
                        recordOffset(input_band_record, metadata.record_length),
                        input_coefficients.data(),
                        input_coefficients.size()
                    )) {
                    std::cerr << "Error reading coefficients for spin " << spin_idx
                              << ", k-point " << kpoint_idx << ", band " << band_idx
                              << std::endl;
                    return -7;
                }

                std::vector<long double> band_norm_diff(fraction_count + 1, 0.0L);
                long double band_norm = 0.0L;

                for (int coeff_idx = 0; coeff_idx < map.original_count; ++coeff_idx) {
                    const std::size_t source_offset =
                        static_cast<std::size_t>(coeff_idx) *
                        coefficient_layout.complex_bytes;
                    const double norm_squared = coefficientNormSquared(
                        input_coefficients.data() + source_offset,
                        coefficient_layout
                    );
                    band_norm += static_cast<long double>(norm_squared);

                    const std::size_t first_fraction =
                        map.first_retained_fraction[static_cast<std::size_t>(coeff_idx)];
                    if (first_fraction < fraction_count) {
                        band_norm_diff[first_fraction] +=
                            static_cast<long double>(norm_squared);
                    }
                }

                total_norm += band_norm;
                long double band_retained = 0.0L;
                for (std::size_t fraction_idx = 0; fraction_idx < fraction_count; ++fraction_idx) {
                    band_retained += band_norm_diff[fraction_idx];
                    retained_norm[fraction_idx] += band_retained;

                    const double band_fraction =
                        band_norm > 0.0L
                            ? static_cast<double>(band_retained / band_norm)
                            : 1.0;
                    min_band_retained_norm[fraction_idx] = std::min(
                        min_band_retained_norm[fraction_idx],
                        band_fraction
                    );
                }
            }
        }
    }

    namespace fs = std::filesystem;
    std::error_code fs_error;
    const std::uint64_t original_file_size = fs::file_size(input_path, fs_error);
    const std::uint64_t record_count =
        2ULL +
        static_cast<std::uint64_t>(metadata.spin_count) *
            static_cast<std::uint64_t>(metadata.kpoint_count) *
            static_cast<std::uint64_t>(metadata.band_count + 1);
    const std::size_t kpoint_header_bytes =
        static_cast<std::size_t>(4 + 3 * metadata.band_count) * sizeof(double);

    analysis.reserve(fraction_count);
    for (std::size_t fraction_idx = 0; fraction_idx < fraction_count; ++fraction_idx) {
        AnalysisPoint point;
        point.cutoff_fraction = fractions[fraction_idx];
        point.cutoff_energy = fractions[fraction_idx] * metadata.cutoff_energy;
        point.total_coefficients = total_coefficients;
        for (const AnalysisPlaneWaveMap& map : maps) {
            point.kept_coefficients +=
                static_cast<std::uint64_t>(map.kept_counts[fraction_idx]) *
                static_cast<std::uint64_t>(metadata.band_count);
        }
        point.coefficient_fraction =
            total_coefficients > 0
                ? static_cast<double>(point.kept_coefficients) /
                    static_cast<double>(total_coefficients)
                : 0.0;
        point.retained_norm_fraction =
            total_norm > 0.0L
                ? static_cast<double>(retained_norm[fraction_idx] / total_norm)
                : 1.0;
        point.lost_norm_fraction = 1.0 - point.retained_norm_fraction;
        point.min_band_retained_norm_fraction = min_band_retained_norm[fraction_idx];

        const std::size_t coefficient_record_bytes =
            static_cast<std::size_t>(max_kept_coefficients[fraction_idx] + 1) *
            coefficient_layout.complex_bytes;
        point.estimated_record_length = static_cast<int>(
            std::max(kpoint_header_bytes, coefficient_record_bytes)
        );
        point.estimated_file_size =
            record_count * static_cast<std::uint64_t>(point.estimated_record_length);
        point.estimated_size_ratio =
            original_file_size > 0
                ? static_cast<double>(point.estimated_file_size) /
                    static_cast<double>(original_file_size)
                : 0.0;

        analysis.push_back(point);
    }

    return 0;
}
