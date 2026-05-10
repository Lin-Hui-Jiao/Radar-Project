#pragma once

#include <optional>
#include <string>
#include <vector>

#include "hiradar/local_tangent_frame.hpp"

namespace hiradar {

struct EnergyGrid2D {
    int width = 0;
    int height = 0;
    double resolution_m = 1.0;
    double fixed_height_m = 0.0;
    LocalTangentFrame origin_frame;
    double x_min_m = 0.0;
    double y_min_m = 0.0;
    std::vector<float> values;

    bool in_bounds(int ix, int iy) const;
    float get(int ix, int iy) const;
    float& at(int ix, int iy);
    size_t point_count() const;
    double max_positive_value(double zero_epsilon = 1e-12) const;
    size_t EstimatedSteadyBytes() const;

    static EnergyGrid2D CreateForTesting(
        int width,
        int height,
        double resolution_m,
        double fixed_height_m,
        const LocalTangentFrame& frame,
        double x_min_m,
        double y_min_m,
        std::vector<float> values
    );

    static EnergyGrid2D LoadFromCsv(
        const std::string& csv_path,
        const std::optional<double>& explicit_resolution_m = std::nullopt
    );
};

}  // namespace hiradar
