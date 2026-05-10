#include "hiradar/energy_grid_2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace hiradar {

namespace {

struct CsvRow {
    double lon = 0.0;
    double lat = 0.0;
    double alt = 0.0;
    float power_density = 0.0f;
};

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        parts.push_back(item);
    }
    return parts;
}

int FindColumnIndex(const std::vector<std::string>& header, const std::string& name) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

double MeanAbsoluteDifference(const std::vector<double>& values) {
    if (values.size() < 2) {
        return 0.0;
    }

    double sum = 0.0;
    size_t count = 0;
    for (size_t i = 1; i < values.size(); ++i) {
        sum += std::abs(values[i] - values[i - 1]);
        ++count;
    }
    return count == 0 ? 0.0 : sum / static_cast<double>(count);
}

double MaxDeviationFromStep(const std::vector<double>& values, double expected_step) {
    if (values.size() < 2) {
        return 0.0;
    }

    double max_deviation = 0.0;
    for (size_t i = 1; i < values.size(); ++i) {
        max_deviation = std::max(max_deviation, std::abs((values[i] - values[i - 1]) - expected_step));
    }
    return max_deviation;
}

}  // namespace

bool EnergyGrid2D::in_bounds(int ix, int iy) const {
    return ix >= 0 && ix < width && iy >= 0 && iy < height;
}

float EnergyGrid2D::get(int ix, int iy) const {
    if (!in_bounds(ix, iy)) {
        throw std::out_of_range("EnergyGrid2D::get index out of range");
    }
    return values[static_cast<size_t>(iy) * width + ix];
}

float& EnergyGrid2D::at(int ix, int iy) {
    if (!in_bounds(ix, iy)) {
        throw std::out_of_range("EnergyGrid2D::at index out of range");
    }
    return values[static_cast<size_t>(iy) * width + ix];
}

size_t EnergyGrid2D::point_count() const {
    return static_cast<size_t>(width) * static_cast<size_t>(height);
}

double EnergyGrid2D::max_positive_value(double zero_epsilon) const {
    double max_value = 0.0;
    for (float value : values) {
        if (value > zero_epsilon) {
            max_value = std::max(max_value, static_cast<double>(value));
        }
    }
    return max_value;
}

size_t EnergyGrid2D::EstimatedSteadyBytes() const {
    return sizeof(EnergyGrid2D) + values.capacity() * sizeof(float);
}

EnergyGrid2D EnergyGrid2D::CreateForTesting(
    int width,
    int height,
    double resolution_m,
    double fixed_height_m,
    const LocalTangentFrame& frame,
    double x_min_m,
    double y_min_m,
    std::vector<float> values
) {
    if (width <= 0 || height <= 0) {
        throw std::runtime_error("CreateForTesting requires positive grid dimensions");
    }
    if (values.size() != static_cast<size_t>(width) * static_cast<size_t>(height)) {
        throw std::runtime_error("CreateForTesting received mismatched value count");
    }

    EnergyGrid2D grid;
    grid.width = width;
    grid.height = height;
    grid.resolution_m = resolution_m;
    grid.fixed_height_m = fixed_height_m;
    grid.origin_frame = frame;
    grid.x_min_m = x_min_m;
    grid.y_min_m = y_min_m;
    grid.values = std::move(values);
    return grid;
}

EnergyGrid2D EnergyGrid2D::LoadFromCsv(
    const std::string& csv_path,
    const std::optional<double>& explicit_resolution_m
) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open CSV file: " + csv_path);
    }

    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("CSV file is empty: " + csv_path);
    }

    const std::vector<std::string> header = SplitCsvLine(header_line);
    const int lon_index = FindColumnIndex(header, "lon");
    const int lat_index = FindColumnIndex(header, "lat");
    const int alt_index = FindColumnIndex(header, "alt");
    const int power_index = FindColumnIndex(header, "power_density");

    if (lon_index < 0 || lat_index < 0 || alt_index < 0 || power_index < 0) {
        throw std::runtime_error("CSV file is missing required columns: lon, lat, alt, power_density");
    }

    std::vector<CsvRow> rows;
    rows.reserve(262144);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        const std::vector<std::string> parts = SplitCsvLine(line);
        const size_t required_columns = static_cast<size_t>(
            std::max(std::max(lon_index, lat_index), std::max(alt_index, power_index)) + 1
        );
        if (parts.size() < required_columns) {
            throw std::runtime_error("CSV row has too few columns: " + line);
        }

        CsvRow row;
        row.lon = std::stod(parts[lon_index]);
        row.lat = std::stod(parts[lat_index]);
        row.alt = std::stod(parts[alt_index]);
        row.power_density = std::stof(parts[power_index]);
        rows.push_back(row);
    }

    if (rows.empty()) {
        throw std::runtime_error("CSV file contains header only: " + csv_path);
    }

    const double alt_tolerance = 1e-4;
    const double fixed_height_m = rows.front().alt;
    double min_lon = rows.front().lon;
    double max_lon = rows.front().lon;
    double min_lat = rows.front().lat;
    double max_lat = rows.front().lat;
    for (const CsvRow& row : rows) {
        if (std::abs(row.alt - fixed_height_m) > alt_tolerance) {
            throw std::runtime_error("CSV contains varying altitudes; expected a fixed-height plane");
        }
        min_lon = std::min(min_lon, row.lon);
        max_lon = std::max(max_lon, row.lon);
        min_lat = std::min(min_lat, row.lat);
        max_lat = std::max(max_lat, row.lat);
    }

    const double row_lat_tolerance = 1e-9;
    std::vector<int> row_offsets;
    std::vector<double> row_lats;
    row_offsets.push_back(0);
    row_lats.push_back(rows.front().lat);

    for (size_t index = 1; index < rows.size(); ++index) {
        if (std::abs(rows[index].lat - row_lats.back()) > row_lat_tolerance) {
            row_offsets.push_back(static_cast<int>(index));
            row_lats.push_back(rows[index].lat);
        }
    }
    row_offsets.push_back(static_cast<int>(rows.size()));

    const int height = static_cast<int>(row_lats.size());
    if (height <= 0) {
        throw std::runtime_error("Failed to recover grid row count from CSV");
    }

    const int width = row_offsets[1] - row_offsets[0];
    if (width <= 0) {
        throw std::runtime_error("Failed to recover grid column count from CSV");
    }

    for (int row = 0; row < height; ++row) {
        const int row_start = row_offsets[row];
        const int row_end = row_offsets[row + 1];
        if (row_end - row_start != width) {
            throw std::runtime_error("CSV row width mismatch detected while recovering grid");
        }
        for (int col = row_start + 1; col < row_end; ++col) {
            if (std::abs(rows[col].lat - rows[row_start].lat) > row_lat_tolerance) {
                throw std::runtime_error("CSV row contains drifting latitudes");
            }
            if (rows[col].lon <= rows[col - 1].lon) {
                throw std::runtime_error("CSV longitudes must be strictly increasing within each row");
            }
        }
        if (row > 0 && row_lats[row] <= row_lats[row - 1]) {
            throw std::runtime_error("CSV latitudes must be strictly increasing across rows");
        }
    }

    std::vector<double> first_row_lons;
    first_row_lons.reserve(width);
    for (int col = row_offsets[0]; col < row_offsets[1]; ++col) {
        first_row_lons.push_back(rows[col].lon);
    }

    std::vector<double> first_column_lats;
    first_column_lats.reserve(height);
    for (int row = 0; row < height; ++row) {
        first_column_lats.push_back(rows[row_offsets[row]].lat);
    }

    const double center_lon = 0.5 * (min_lon + max_lon);
    const double center_lat = 0.5 * (min_lat + max_lat);
    LocalTangentFrame frame(center_lon, center_lat, fixed_height_m);

    std::vector<double> east_coords;
    east_coords.reserve(width);
    for (double lon : first_row_lons) {
        east_coords.push_back(frame.ToENU(lon, first_column_lats.front(), fixed_height_m).x);
    }

    std::vector<double> north_coords;
    north_coords.reserve(height);
    for (double lat : first_column_lats) {
        north_coords.push_back(frame.ToENU(first_row_lons.front(), lat, fixed_height_m).y);
    }

    const double east_step = MeanAbsoluteDifference(east_coords);
    const double north_step = MeanAbsoluteDifference(north_coords);
    double inferred_resolution = 0.5 * (east_step + north_step);
    if (explicit_resolution_m.has_value()) {
        inferred_resolution = *explicit_resolution_m;
    }

    if (inferred_resolution <= 0.0) {
        throw std::runtime_error("Inferred grid resolution is invalid");
    }

    const double step_tolerance = std::max(1e-3, 0.05 * inferred_resolution);
    if (std::abs(east_step - north_step) > step_tolerance) {
        throw std::runtime_error("Recovered east/north resolution mismatch exceeds tolerance");
    }
    if (MaxDeviationFromStep(east_coords, east_step) > step_tolerance) {
        throw std::runtime_error("Recovered longitude spacing drift exceeds tolerance");
    }
    if (MaxDeviationFromStep(north_coords, north_step) > step_tolerance) {
        throw std::runtime_error("Recovered latitude spacing drift exceeds tolerance");
    }
    if (explicit_resolution_m.has_value()) {
        if (std::abs(east_step - *explicit_resolution_m) > step_tolerance ||
            std::abs(north_step - *explicit_resolution_m) > step_tolerance) {
            throw std::runtime_error("Explicit resolution does not match CSV-derived spacing");
        }
    }

    EnergyGrid2D grid;
    grid.width = width;
    grid.height = height;
    grid.resolution_m = inferred_resolution;
    grid.fixed_height_m = fixed_height_m;
    grid.origin_frame = frame;
    grid.x_min_m = east_coords.front();
    grid.y_min_m = north_coords.front();
    grid.values.assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0.0f);

    for (int row = 0; row < height; ++row) {
        const int row_start = row_offsets[row];
        for (int col = 0; col < width; ++col) {
            const CsvRow& sample = rows[row_start + col];
            grid.values[static_cast<size_t>(row) * grid.width + col] = sample.power_density;
        }
    }

    return grid;
}

}  // namespace hiradar
