#include "reader.hpp"

int readMetadata(WavecarMetadata& meta_data, std::ifstream& file) {
    static_assert(sizeof(double) == 8, "Size of double must be 8 bytes");
    // read the metadata
    double buffer;
    file.read(reinterpret_cast<char*>(&buffer), 8); // 8 = size of double
    meta_data.record_length = static_cast<int>(buffer);
    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.spin_count = static_cast<int>(buffer);
    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.RTAG_TAG = static_cast<int>(buffer);

    file.seekg(meta_data.record_length, std::ios::beg); // move to header part

    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.kpoint_count = static_cast<int>(buffer);
    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.band_count = static_cast<int>(buffer);
    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.cutoff_energy = buffer;
    
    // read lattice vectors
    Vec3 a, b, c;
    file.read(reinterpret_cast<char*>(&a.x), 8);
    file.read(reinterpret_cast<char*>(&a.y), 8);
    file.read(reinterpret_cast<char*>(&a.z), 8);
    file.read(reinterpret_cast<char*>(&b.x), 8);
    file.read(reinterpret_cast<char*>(&b.y), 8);
    file.read(reinterpret_cast<char*>(&b.z), 8);
    file.read(reinterpret_cast<char*>(&c.x), 8);
    file.read(reinterpret_cast<char*>(&c.y), 8);
    file.read(reinterpret_cast<char*>(&c.z), 8);

    meta_data.lattice = Cell(a, b, c);

    file.read(reinterpret_cast<char*>(&buffer), 8);
    meta_data.fermi_energy = buffer;

    // check if the read was successful
    if (!file) {
        return -1;
    }
    return 0;
}

// TODO: Handle complex eigenvalues if needed
int readEigenvalues(std::unique_ptr<double[]>& eigenvalues, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file) {
    //int kpoint_data_size = (meta_data.band_count + 1) * meta_data.record_length; // nband coeffs + 1 eig/occ
    //int spin_data_size = kpoint_data_size * meta_data.kpoint_count;
    //int starting_byte_offset = spin_idx * spin_data_size + kpoint_idx * kpoint_data_size + 2 * meta_data.record_length + 4 * sizeof(double); // + offset for header & meta + coeffcount & kpoint
    int starting_byte_offset = (spin_idx * meta_data.kpoint_count + kpoint_idx) * (meta_data.band_count + 1) * meta_data.record_length + 2 * meta_data.record_length;
    starting_byte_offset += 4 * sizeof(double); // skip number of coeffs and kpoint
    file.seekg(starting_byte_offset, std::ios::beg);
    if(!file) {
        return -1; // error seeking in file
    }
    for (int i = 0; i < meta_data.band_count; ++i)
    {
        double real, imag;
        file.read(reinterpret_cast<char*>(&real), sizeof(double));
        file.read(reinterpret_cast<char*>(&imag), sizeof(double));
        if (imag != 0.0) {
            // Handle complex eigenvalues if needed
            std::cerr << "Complex eigenvalues are not supported in this implementation." << std::endl;
            return -2; // complex eigenvalues not supported
        }
        eigenvalues[i] = real;
        file.seekg(sizeof(double), std::ios::cur); // skip occupancy
    }
    if (!file) {
        return -3; // error reading eigenvalues
    }
    return 0; // success
}

int getBandCoefficientsSize(int& num_coeffs, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file) {
    int kpoint_data_size = (meta_data.band_count + 1) * meta_data.record_length; // nband coeffs + 1 eig/occ
    int spin_data_size = kpoint_data_size * meta_data.kpoint_count;
    int starting_byte_offset = spin_idx * spin_data_size + kpoint_idx * kpoint_data_size + 2 * meta_data.record_length;
    double coeff_count;
    file.seekg(starting_byte_offset, std::ios::beg);
    if(!file) {
        return -1; // error seeking in file
    }
    file.read(reinterpret_cast<char*>(&coeff_count), sizeof(double));
    if (!file) {
        return -2; // error reading coefficient count
    }
    num_coeffs = static_cast<int>(coeff_count);
    if (num_coeffs <= 0) {
        return -3; // invalid number of coefficients
    }
    return 0;
}

int readBandCoefficients(std::vector<double>& coeffs, int band_idx, int kpoint_idx, int spin_idx, const WavecarMetadata& meta_data, std::ifstream& file){
    int kpoint_data_size = (meta_data.band_count + 1) * meta_data.record_length; // nband coeffs + 1 eig/occ
    int spin_data_size = kpoint_data_size * meta_data.kpoint_count;
    int starting_byte_offset = spin_idx * spin_data_size + kpoint_idx * kpoint_data_size + 2 * meta_data.record_length;
    
    file.seekg(starting_byte_offset, std::ios::beg);
    if(!file) {
        return -1; // error seeking in file
    }
    double num_coeffs;
    file.read(reinterpret_cast<char*>(&num_coeffs), sizeof(double));
    int band_coeffs_offset = starting_byte_offset + band_idx * meta_data.record_length;
    file.seekg(band_coeffs_offset, std::ios::beg);
    // read coefficients
    file.read(reinterpret_cast<char*>(coeffs.data()), num_coeffs * sizeof(double));
    if (!file) {
        return -2; // error reading coefficients
    }
    return 0;
}