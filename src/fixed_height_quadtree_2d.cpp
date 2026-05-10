// This file is the CLI entrypoint for the 2D fixed-height adaptive quadtree prototype.
// It does not implement the quadtree algorithm itself. Instead, it is responsible for
// orchestrating the full workflow:
// 1. parse command-line arguments
// 2. load or build a regular fixed-height energy grid
// 3. configure split thresholds
// 4. build the adaptive quadtree
// 5. export experiment outputs
// 6. optionally run a few point queries for quick validation

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>

#include "hiradar/adaptive_quadtree_2d.hpp"
#include "hiradar/fixed_height_field_builder_2d.hpp"

namespace {

using hiradar::AdaptiveQuadtree2D;
using hiradar::AdaptiveSplitConfig;
using hiradar::EnergyGrid2D;
using hiradar::FixedHeightFieldBuilder2D;
using hiradar::MeshFieldBuildConfig;

// ParsedArgs keeps command-line input split into two categories:
// - options: arguments that look like "--key value"
// - flags: arguments that look like "--flag" with no following value
// This lightweight representation is enough for the current prototype.
struct ParsedArgs {
    std::map<std::string, std::string> options;
    std::set<std::string> flags;
};

// Minimal argument parser.
// Every token must start with "--".
// If the next token does not start with "--", it is treated as the value of the current option.
// Otherwise, the current token is treated as a standalone flag.
ParsedArgs ParseArgs(int argc, char** argv) {
    ParsedArgs args;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key.rfind("--", 0) != 0) {
            throw std::runtime_error("Unexpected positional argument: " + key);
        }
        if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
            args.options[key] = argv[++i];
        } else {
            args.flags.insert(key);
        }
    }
    return args;
}

// Small helper functions used in main() to keep the control flow readable.
bool HasFlag(const ParsedArgs& args, const std::string& flag) {
    return args.flags.find(flag) != args.flags.end();
}

bool HasOption(const ParsedArgs& args, const std::string& key) {
    return args.options.find(key) != args.options.end();
}

std::string GetString(const ParsedArgs& args, const std::string& key, const std::string& fallback) {
    auto it = args.options.find(key);
    return it == args.options.end() ? fallback : it->second;
}

double GetDouble(const ParsedArgs& args, const std::string& key, double fallback) {
    auto it = args.options.find(key);
    return it == args.options.end() ? fallback : std::stod(it->second);
}

int GetInt(const ParsedArgs& args, const std::string& key, int fallback) {
    auto it = args.options.find(key);
    return it == args.options.end() ? fallback : std::stoi(it->second);
}

// Prints supported invocation forms.
// The executable currently supports exactly two input modes:
// - CSV mode: consume an already-generated regular truth field
// - mesh mode: compute the regular truth field from the RTree-based mesh path first
void PrintUsage() {
    std::cout
        << "Usage:\n"
        << "  fixed_height_quadtree_2d --input-csv <path> [options]\n"
        << "  fixed_height_quadtree_2d --build-mesh-field [options]\n\n"
        << "Core options:\n"
        << "  --output-prefix <prefix>           Output prefix for *_leaves.csv and *_summary.txt\n"
        << "  --gradient-threshold <value>       Relative corner-gradient split threshold\n"
        << "  --error-threshold <value>          Absolute center-error split threshold\n"
        << "  --min-block-size <value>           Minimum quadtree block size in cells\n"
        << "  --query-ix <int> --query-iy <int>  Query one grid index after building\n"
        << "  --query-x <m> --query-y <m>        Query one local-meter location after building\n"
        << "  --query-lon <deg> --query-lat <deg> Query one geodetic location after building\n\n"
        << "Mesh-field options:\n"
        << "  --rtree <path>\n"
        << "  --lon-min <deg> --lon-max <deg>\n"
        << "  --lat-min <deg> --lat-max <deg>\n"
        << "  --fixed-height <m>\n"
        << "  --grid-resolution <m>\n"
        << "  --radar-lon <deg> --radar-lat <deg>\n"
        << "  --radar-relative-alt <m>\n";
}

// Writes a compact experiment summary for later inspection and reproducibility.
// The summary records:
// - original and padded grid dimensions
// - leaf statistics
// - the thresholds actually used in this run
// - the local ENU frame origin and the grid origin in local meters
void WriteSummaryFile(
    const std::string& output_path,
    const EnergyGrid2D& grid,
    const AdaptiveSplitConfig& split_config,
    const AdaptiveQuadtree2D& quadtree
) {
    const hiradar::AdaptiveQuadtreeSummary& summary = quadtree.summary();
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open summary output: " + output_path);
    }

    file << std::fixed << std::setprecision(6);
    file << "original_width=" << summary.original_width << "\n";
    file << "original_height=" << summary.original_height << "\n";
    file << "padded_size=" << summary.padded_size << "\n";
    file << "original_points=" << summary.original_points << "\n";
    file << "valid_points=" << summary.valid_points << "\n";
    file << "leaf_count=" << summary.leaf_count << "\n";
    file << "zero_leaf_count=" << summary.zero_leaf_count << "\n";
    file << "smooth_leaf_count=" << summary.smooth_leaf_count << "\n";
    file << "boundary_leaf_count=" << summary.boundary_leaf_count << "\n";
    file << "grid_resolution_m=" << grid.resolution_m << "\n";
    file << "fixed_height_m=" << grid.fixed_height_m << "\n";
    file << "global_max_positive=" << summary.global_max_positive << "\n";
    file << "gradient_threshold_rel=" << split_config.gradient_threshold_rel << "\n";
    file << "center_error_threshold_abs=" << split_config.center_error_threshold_abs << "\n";
    file << "min_block_size=" << split_config.min_block_size << "\n";
    file << "zero_epsilon=" << split_config.zero_epsilon << "\n";
    file << "origin_lon_deg=" << grid.origin_frame.origin_lon_deg() << "\n";
    file << "origin_lat_deg=" << grid.origin_frame.origin_lat_deg() << "\n";
    file << "origin_h_m=" << grid.origin_frame.origin_h_m() << "\n";
    file << "x_min_m=" << grid.x_min_m << "\n";
    file << "y_min_m=" << grid.y_min_m << "\n";
}

// Prints the same key information to stdout for immediate feedback after a run.
// This is useful during iterative experiments, because you can quickly inspect:
// - recovered grid size
// - whether padding was needed
// - how many leaves were created
// - how aggressive the split parameters were
void PrintSummaryToStdout(
    const EnergyGrid2D& grid,
    const AdaptiveSplitConfig& split_config,
    const AdaptiveQuadtree2D& quadtree
) {
    const hiradar::AdaptiveQuadtreeSummary& summary = quadtree.summary();
    std::cout << "Grid: " << summary.original_width << " x " << summary.original_height
              << " (padded to " << summary.padded_size << ")\n";
    std::cout << "Resolution: " << grid.resolution_m << " m, fixed height: " << grid.fixed_height_m << " m\n";
    std::cout << "Leaves: " << summary.leaf_count
              << " [zero=" << summary.zero_leaf_count
              << ", smooth=" << summary.smooth_leaf_count
              << ", boundary=" << summary.boundary_leaf_count << "]\n";
    std::cout << "Thresholds: gradient=" << split_config.gradient_threshold_rel
              << ", center_error=" << split_config.center_error_threshold_abs
              << ", min_block=" << split_config.min_block_size << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        // Step 1: parse command-line arguments.
        // This keeps the rest of main() focused on the actual data-processing pipeline.
        const ParsedArgs args = ParseArgs(argc, argv);
        if (argc == 1 || HasFlag(args, "--help")) {
            PrintUsage();
            return 0;
        }

        // Step 2: determine which input mode is active.
        // Exactly one mode must be selected:
        // - CSV mode means the regular truth field already exists
        // - mesh mode means we first generate the regular truth field from RTree + radar
        const bool csv_mode = HasOption(args, "--input-csv");
        const bool mesh_mode = HasFlag(args, "--build-mesh-field");
        if (csv_mode == mesh_mode) {
            throw std::runtime_error("Choose exactly one input mode: --input-csv or --build-mesh-field");
        }

        EnergyGrid2D grid;
        if (csv_mode) {
            // CSV-backed path:
            // Load a precomputed fixed-height regular grid from an existing CSV file.
            // This is the preferred path when the truth field is already available.
            const std::optional<double> explicit_resolution = HasOption(args, "--grid-resolution")
                ? std::optional<double>(GetDouble(args, "--grid-resolution", 1.0))
                : std::nullopt;
            grid = EnergyGrid2D::LoadFromCsv(GetString(args, "--input-csv", ""), explicit_resolution);
        } else {
            // Mesh-backed path:
            // Build a regular fixed-height truth field on the fly using the mesh-index-based
            // occlusion path and the radar power-density model.
            MeshFieldBuildConfig config;
            config.rtree_file = GetString(args, "--rtree", config.rtree_file);
            config.lon_min = GetDouble(args, "--lon-min", config.lon_min);
            config.lon_max = GetDouble(args, "--lon-max", config.lon_max);
            config.lat_min = GetDouble(args, "--lat-min", config.lat_min);
            config.lat_max = GetDouble(args, "--lat-max", config.lat_max);
            config.fixed_height_m = GetDouble(args, "--fixed-height", config.fixed_height_m);
            config.resolution_m = GetDouble(args, "--grid-resolution", config.resolution_m);
            config.radar_lon_deg = GetDouble(args, "--radar-lon", config.radar_lon_deg);
            config.radar_lat_deg = GetDouble(args, "--radar-lat", config.radar_lat_deg);
            config.radar_relative_alt_m = GetDouble(args, "--radar-relative-alt", config.radar_relative_alt_m);
            grid = FixedHeightFieldBuilder2D::BuildFromMesh(config);
        }

        // Step 3: configure quadtree split parameters.
        // If the user does not provide an absolute center-error threshold, we derive one
        // from the current grid's maximum positive energy value so that the default stays
        // roughly proportional to the signal scale.
        AdaptiveSplitConfig split_config;
        split_config.gradient_threshold_rel = GetDouble(args, "--gradient-threshold", split_config.gradient_threshold_rel);
        split_config.min_block_size = GetInt(args, "--min-block-size", split_config.min_block_size); // 最小能细分到多少块
        if (HasOption(args, "--error-threshold")) {
            split_config.center_error_threshold_abs = GetDouble(
                args,
                "--error-threshold",
                split_config.center_error_threshold_abs
            );
        } else {
            // 如果用户没传，就自动设置成当前网格最大正能量值的 5%
            split_config.center_error_threshold_abs = 0.05 * grid.max_positive_value(split_config.zero_epsilon);
        }

        // Step 4: build the adaptive quadtree.
        // Internally this stage performs padding to the next power-of-two square domain,
        // computes block statistics, recursively splits where needed, and stores the final leaves.
        AdaptiveQuadtree2D quadtree;
        quadtree.Build(grid, split_config);

        // Step 5: export experiment artifacts.
        // - *_leaves.csv stores the final adaptive blocks
        // - *_summary.txt stores build metadata and thresholds
        const std::string output_prefix = GetString(args, "--output-prefix", "fixed_height_quadtree_2d");
        quadtree.ExportLeavesCsv(output_prefix + "_leaves.csv");
        WriteSummaryFile(output_prefix + "_summary.txt", grid, split_config, quadtree);
        PrintSummaryToStdout(grid, split_config, quadtree);

        // Step 6: optional point queries.
        // These are mainly for smoke tests and manual debugging after a run.
        // The three query entrances correspond to three coordinate views:
        // - grid indices (ix, iy)
        // - local ENU meters (x, y)
        // - geodetic coordinates (lon, lat)
        if (HasOption(args, "--query-ix") && HasOption(args, "--query-iy")) {
            const auto value = quadtree.QueryByIndex(
                GetInt(args, "--query-ix", 0),
                GetInt(args, "--query-iy", 0)
            );
            std::cout << "QueryByIndex: " << (value.has_value() ? std::to_string(*value) : "miss") << "\n";
        }
        if (HasOption(args, "--query-x") && HasOption(args, "--query-y")) {
            const auto value = quadtree.QueryByLocalMeters(
                GetDouble(args, "--query-x", 0.0),
                GetDouble(args, "--query-y", 0.0)
            );
            std::cout << "QueryByLocalMeters: " << (value.has_value() ? std::to_string(*value) : "miss") << "\n";
        }
        if (HasOption(args, "--query-lon") && HasOption(args, "--query-lat")) {
            const auto value = quadtree.QueryByGeodetic(
                GetDouble(args, "--query-lon", 0.0),
                GetDouble(args, "--query-lat", 0.0)
            );
            std::cout << "QueryByGeodetic: " << (value.has_value() ? std::to_string(*value) : "miss") << "\n";
        }

        return 0;
    } catch (const std::exception& ex) {
        // Funnel all operational errors to one place so the CLI exits cleanly with a readable message.
        std::cerr << "fixed_height_quadtree_2d error: " << ex.what() << std::endl;
        return 1;
    }
}
