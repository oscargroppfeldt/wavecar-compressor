#pragma once
#ifndef READER_HPP
#define READER_HPP

#include "util_structs.hpp"
#include <fstream>
#include <iostream>
#include <memory>

int readMetadata(WavecarMetadata& meta_data, std::ifstream& file);

int readEigenvalues(std::unique_ptr<double[]>& eigenvalues, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file);

int getBandCoefficientsSize(int& num_coeffs, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file);

int readBandCoefficients(std::vector<double>& coeffs, int band_idx, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file);

#endif // READER_HPP