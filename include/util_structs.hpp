#pragma once
#ifndef UTIL_STRUCTS_HPP
#define UTIL_STRUCTS_HPP

#include <cmath>
#include <vector>

struct Vec3{
    double x;
    double y;
    double z;

    Vec3(double x_val = 0.0, double y_val = 0.0, double z_val = 0.0)
        : x(x_val), y(y_val), z(z_val) {}

    double dot(Vec3 other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    Vec3 cross(Vec3 other) const {
        return Vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
};

struct Cell{
    Vec3 a;
    Vec3 b;
    Vec3 c;

    double vol;

    Cell(Vec3 a_val = Vec3(), Vec3 b_val = Vec3(), Vec3 c_val = Vec3())
        : a(a_val), b(b_val), c(c_val), vol(0.0) {vol = volume();}

    double volume() {
        if (vol != 0.0) {
            return vol; // Return cached volume if already calculated
        }
        vol = a.dot(b.cross(c));
        return vol;
    }
};

struct WavecarMetadata {
    int record_length; // in bytes
    int spin_count;
    int RTAG_TAG;
    int kpoint_count;
    int band_count;
    double cutoff_energy;
    Cell lattice;
    double fermi_energy;
};


#endif // UTIL_STRUCTS_HPP
