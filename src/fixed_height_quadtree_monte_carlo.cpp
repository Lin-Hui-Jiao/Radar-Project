#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "hiradar/adaptive_quadtree_2d.hpp"

namespace {

using hiradar::AdaptiveQuadtree2D;
using hiradar::AdaptiveQuadtreeMemoryFootprint;
using hiradar::AdaptiveQuadtreeSummary;
using hiradar::AdaptiveSplitConfig;
using hiradar::EnergyGrid2D;

struct ParsedArgs {
    std::map<std::string, std::string> options;
    std::set<std::string> flags;
};

struct IndexRegion {
    int ix_min = 0;
    int ix_max = 0;
    int iy_min = 0;
    int iy_max = 0;
    bool valid = true;
    std::string invalid_reason;
};

enum class SamplingMode {
    Cell,
    Continuous,
};

enum class DistributionMode {
    Uniform,
    Boundary,
};

struct QuerySample {
    int ix = 0;
    int iy = 0;
    double x_m = 0.0;
    double y_m = 0.0;
};

struct ScenarioDefinition {
    std::string name;
    SamplingMode sampling_mode = SamplingMode::Cell;
    DistributionMode distribution = DistributionMode::Uniform;
    size_t requested_samples = 0;
};

struct AccuracyStats {
    size_t sample_count = 0;
    size_t hit_count_baseline = 0;
    size_t hit_count_tree = 0;
    double mae = 0.0;
    double rmse = 0.0;
    double max_abs_error = 0.0;
    double p95_abs_error = 0.0;
};

struct TimingStats {
    std::vector<double> baseline_round_ns_per_query;
    std::vector<double> tree_round_ns_per_query;
    double baseline_mean_ns_per_query = 0.0;
    double baseline_median_ns_per_query = 0.0;
    double tree_mean_ns_per_query = 0.0;
    double tree_median_ns_per_query = 0.0;
    double speedup_ratio = 0.0;
};

struct ScenarioMetrics {
    std::string scenario_name;
    std::string sampling_mode;
    std::string distribution;
    size_t sample_count = 0;
    uint64_t seed = 0;
    size_t hit_count_baseline = 0;
    size_t hit_count_tree = 0;
    double mae = 0.0;
    double rmse = 0.0;
    double max_abs_error = 0.0;
    double p95_abs_error = 0.0;
    double baseline_mean_ns_per_query = 0.0;
    double baseline_median_ns_per_query = 0.0;
    double tree_mean_ns_per_query = 0.0;
    double tree_median_ns_per_query = 0.0;
    double speedup_ratio = 0.0;
    size_t grid_steady_bytes = 0;
    size_t tree_steady_bytes = 0;
    size_t tree_peak_incremental_bytes = 0;
    size_t tree_total_peak_bytes = 0;
    double tree_build_ms = 0.0;
    size_t leaf_count = 0;
    size_t zero_leaf_count = 0;
    size_t smooth_leaf_count = 0;
    size_t boundary_leaf_count = 0;
    bool skipped = false;
    std::string skip_reason;
};

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

bool HasOption(const ParsedArgs& args, const std::string& key) {
    return args.options.find(key) != args.options.end();
}

bool HasFlag(const ParsedArgs& args, const std::string& key) {
    return args.flags.find(key) != args.flags.end();
}

std::string GetString(const ParsedArgs& args, const std::string& key, const std::string& fallback) {
    const auto it = args.options.find(key);
    return it == args.options.end() ? fallback : it->second;
}

double GetDouble(const ParsedArgs& args, const std::string& key, double fallback) {
    const auto it = args.options.find(key);
    return it == args.options.end() ? fallback : std::stod(it->second);
}

int GetInt(const ParsedArgs& args, const std::string& key, int fallback) {
    const auto it = args.options.find(key);
    return it == args.options.end() ? fallback : std::stoi(it->second);
}

uint64_t GetUInt64(const ParsedArgs& args, const std::string& key, uint64_t fallback) {
    const auto it = args.options.find(key);
    return it == args.options.end() ? fallback : static_cast<uint64_t>(std::stoull(it->second));
}

void PrintUsage() {
    std::cout
        << "Usage:\n"
        << "  fixed_height_quadtree_monte_carlo --input-csv <path> [options]\n\n"
        << "Core options:\n"
        << "  --output-prefix <prefix>\n"
        << "  --gradient-threshold <value>\n"
        << "  --error-threshold <value>\n"
        << "  --min-block-size <value>\n"
        << "  --full-samples <int>\n"
        << "  --boundary-samples <int>\n"
        << "  --warmup-rounds <int>\n"
        << "  --timing-rounds <int>\n"
        << "  --seed <uint64>\n"
        << "  --ix-min <int> --ix-max <int> --iy-min <int> --iy-max <int>\n";
}

std::string SamplingModeToString(SamplingMode mode) {
    return mode == SamplingMode::Cell ? "cell" : "continuous";
}

std::string DistributionModeToString(DistributionMode mode) {
    return mode == DistributionMode::Uniform ? "uniform" : "boundary";
}

double UniformUnit(std::mt19937_64& rng) {
    return std::generate_canonical<double, 53>(rng);
}

double Mean(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / static_cast<double>(values.size());
}

double Median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if ((values.size() % 2) == 0) {
        return 0.5 * (values[mid - 1] + values[mid]);
    }
    return values[mid];
}

double Percentile95(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t index = static_cast<size_t>(std::ceil(0.95 * static_cast<double>(values.size()))) - 1;
    return values[std::min(index, values.size() - 1)];
}

IndexRegion ResolveRegion(const ParsedArgs& args, const EnergyGrid2D& grid) {
    const bool any_region_option =
        HasOption(args, "--ix-min") ||
        HasOption(args, "--ix-max") ||
        HasOption(args, "--iy-min") ||
        HasOption(args, "--iy-max");

    IndexRegion region;
    if (!any_region_option) {
        region.ix_min = 0;
        region.ix_max = grid.width;
        region.iy_min = 0;
        region.iy_max = grid.height;
        return region;
    }

    if (!HasOption(args, "--ix-min") || !HasOption(args, "--ix-max") ||
        !HasOption(args, "--iy-min") || !HasOption(args, "--iy-max")) {
        region.valid = false;
        region.invalid_reason = "partial_region_arguments";
        return region;
    }

    region.ix_min = GetInt(args, "--ix-min", 0);
    region.ix_max = GetInt(args, "--ix-max", 0);
    region.iy_min = GetInt(args, "--iy-min", 0);
    region.iy_max = GetInt(args, "--iy-max", 0);

    if (region.ix_min < 0 || region.iy_min < 0 ||
        region.ix_max > grid.width || region.iy_max > grid.height ||
        region.ix_min >= region.ix_max || region.iy_min >= region.iy_max) {
        region.valid = false;
        region.invalid_reason = "invalid_region";
    }
    return region;
}

std::vector<std::pair<int, int>> BuildBoundaryCells(
    const EnergyGrid2D& grid,
    const IndexRegion& region,
    const AdaptiveSplitConfig& split_config
) {
    std::vector<std::pair<int, int>> cells;
    for (int iy = region.iy_min; iy < region.iy_max; ++iy) {
        for (int ix = region.ix_min; ix < region.ix_max; ++ix) {
            int zero_count = 0;
            int positive_count = 0;
            float max_value = grid.get(ix, iy);
            float min_value = grid.get(ix, iy);
            for (int ny = std::max(0, iy - 1); ny <= std::min(grid.height - 1, iy + 1); ++ny) {
                for (int nx = std::max(0, ix - 1); nx <= std::min(grid.width - 1, ix + 1); ++nx) {
                    const float value = grid.get(nx, ny);
                    max_value = std::max(max_value, value);
                    min_value = std::min(min_value, value);
                    if (value <= split_config.zero_epsilon) {
                        ++zero_count;
                    } else {
                        ++positive_count;
                    }
                }
            }

            const bool zero_nonzero_mixed = zero_count > 0 && positive_count > 0;
            const double rel_range = (max_value - min_value) /
                std::max(static_cast<double>(max_value), split_config.zero_epsilon);
            if (zero_nonzero_mixed || rel_range > split_config.gradient_threshold_rel) {
                cells.emplace_back(ix, iy);
            }
        }
    }
    return cells;
}

std::optional<float> BaselineQueryByIndex(const EnergyGrid2D& grid, int ix, int iy) {
    if (!grid.in_bounds(ix, iy)) {
        return std::nullopt;
    }
    return grid.get(ix, iy);
}

std::optional<float> BaselineQueryByLocalMeters(const EnergyGrid2D& grid, double x_m, double y_m) {
    if (grid.resolution_m <= 0.0) {
        return std::nullopt;
    }
    const int ix = static_cast<int>(std::floor((x_m - grid.x_min_m) / grid.resolution_m));
    const int iy = static_cast<int>(std::floor((y_m - grid.y_min_m) / grid.resolution_m));
    return BaselineQueryByIndex(grid, ix, iy);
}

std::vector<QuerySample> GenerateSamples(
    const ScenarioDefinition& scenario,
    const EnergyGrid2D& grid,
    const IndexRegion& region,
    const std::vector<std::pair<int, int>>& boundary_cells,
    std::mt19937_64& rng
) {
    std::vector<QuerySample> samples;
    samples.reserve(scenario.requested_samples);

    if (scenario.distribution == DistributionMode::Boundary && boundary_cells.empty()) {
        return samples;
    }

    const double resolution = grid.resolution_m;
    const double x_base = grid.x_min_m;
    const double y_base = grid.y_min_m;

    if (scenario.distribution == DistributionMode::Uniform && scenario.sampling_mode == SamplingMode::Cell) {
        std::uniform_int_distribution<int> dist_x(region.ix_min, region.ix_max - 1);
        std::uniform_int_distribution<int> dist_y(region.iy_min, region.iy_max - 1);
        for (size_t i = 0; i < scenario.requested_samples; ++i) {
            QuerySample sample;
            sample.ix = dist_x(rng);
            sample.iy = dist_y(rng);
            sample.x_m = x_base + (static_cast<double>(sample.ix) + 0.5) * resolution;
            sample.y_m = y_base + (static_cast<double>(sample.iy) + 0.5) * resolution;
            samples.push_back(sample);
        }
        return samples;
    }

    if (scenario.distribution == DistributionMode::Uniform && scenario.sampling_mode == SamplingMode::Continuous) {
        const double x_min = x_base + static_cast<double>(region.ix_min) * resolution;
        const double x_max = x_base + static_cast<double>(region.ix_max) * resolution;
        const double y_min = y_base + static_cast<double>(region.iy_min) * resolution;
        const double y_max = y_base + static_cast<double>(region.iy_max) * resolution;
        for (size_t i = 0; i < scenario.requested_samples; ++i) {
            QuerySample sample;
            sample.x_m = x_min + UniformUnit(rng) * (x_max - x_min);
            sample.y_m = y_min + UniformUnit(rng) * (y_max - y_min);
            sample.ix = static_cast<int>(std::floor((sample.x_m - x_base) / resolution));
            sample.iy = static_cast<int>(std::floor((sample.y_m - y_base) / resolution));
            samples.push_back(sample);
        }
        return samples;
    }

    std::uniform_int_distribution<size_t> boundary_pick(0, boundary_cells.size() - 1);
    for (size_t i = 0; i < scenario.requested_samples; ++i) {
        const auto& cell = boundary_cells[boundary_pick(rng)];
        QuerySample sample;
        sample.ix = cell.first;
        sample.iy = cell.second;
        if (scenario.sampling_mode == SamplingMode::Cell) {
            sample.x_m = x_base + (static_cast<double>(sample.ix) + 0.5) * resolution;
            sample.y_m = y_base + (static_cast<double>(sample.iy) + 0.5) * resolution;
        } else {
            sample.x_m = x_base + static_cast<double>(sample.ix) * resolution + UniformUnit(rng) * resolution;
            sample.y_m = y_base + static_cast<double>(sample.iy) * resolution + UniformUnit(rng) * resolution;
        }
        samples.push_back(sample);
    }
    return samples;
}

AccuracyStats EvaluateAccuracy(
    const ScenarioDefinition& scenario,
    const EnergyGrid2D& grid,
    const AdaptiveQuadtree2D& tree,
    const std::vector<QuerySample>& samples
) {
    AccuracyStats stats;
    stats.sample_count = samples.size();
    std::vector<double> abs_errors;
    abs_errors.reserve(samples.size());

    double sum_abs = 0.0;
    double sum_sq = 0.0;
    for (const QuerySample& sample : samples) {
        const auto baseline = scenario.sampling_mode == SamplingMode::Cell
            ? BaselineQueryByIndex(grid, sample.ix, sample.iy)
            : BaselineQueryByLocalMeters(grid, sample.x_m, sample.y_m);
        const auto tree_value = scenario.sampling_mode == SamplingMode::Cell
            ? tree.QueryByIndex(sample.ix, sample.iy)
            : tree.QueryByLocalMeters(sample.x_m, sample.y_m);

        if (!baseline.has_value() || !tree_value.has_value()) {
            throw std::runtime_error("Benchmark query miss encountered inside a valid sample set");
        }

        ++stats.hit_count_baseline;
        ++stats.hit_count_tree;

        const double abs_error = std::abs(static_cast<double>(*tree_value) - static_cast<double>(*baseline));
        abs_errors.push_back(abs_error);
        sum_abs += abs_error;
        sum_sq += abs_error * abs_error;
        stats.max_abs_error = std::max(stats.max_abs_error, abs_error);
    }

    if (!samples.empty()) {
        const double count = static_cast<double>(samples.size());
        stats.mae = sum_abs / count;
        stats.rmse = std::sqrt(sum_sq / count);
        stats.p95_abs_error = Percentile95(std::move(abs_errors));
    }
    return stats;
}

double RunBaselineRound(
    const ScenarioDefinition& scenario,
    const EnergyGrid2D& grid,
    const std::vector<QuerySample>& samples
) {
    volatile double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (const QuerySample& sample : samples) {
        const auto value = scenario.sampling_mode == SamplingMode::Cell
            ? BaselineQueryByIndex(grid, sample.ix, sample.iy)
            : BaselineQueryByLocalMeters(grid, sample.x_m, sample.y_m);
        if (!value.has_value()) {
            throw std::runtime_error("Baseline benchmark query missed during timing run");
        }
        checksum += *value;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    return samples.empty() ? 0.0 : static_cast<double>(elapsed_ns) / static_cast<double>(samples.size());
}

double RunTreeRound(
    const ScenarioDefinition& scenario,
    const AdaptiveQuadtree2D& tree,
    const std::vector<QuerySample>& samples
) {
    volatile double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (const QuerySample& sample : samples) {
        const auto value = scenario.sampling_mode == SamplingMode::Cell
            ? tree.QueryByIndex(sample.ix, sample.iy)
            : tree.QueryByLocalMeters(sample.x_m, sample.y_m);
        if (!value.has_value()) {
            throw std::runtime_error("Tree benchmark query missed during timing run");
        }
        checksum += *value;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    return samples.empty() ? 0.0 : static_cast<double>(elapsed_ns) / static_cast<double>(samples.size());
}

TimingStats EvaluateTiming(
    const ScenarioDefinition& scenario,
    const EnergyGrid2D& grid,
    const AdaptiveQuadtree2D& tree,
    const std::vector<QuerySample>& samples,
    int warmup_rounds,
    int timing_rounds
) {
    TimingStats stats;
    for (int i = 0; i < warmup_rounds; ++i) {
        (void)RunBaselineRound(scenario, grid, samples);
    }
    for (int i = 0; i < warmup_rounds; ++i) {
        (void)RunTreeRound(scenario, tree, samples);
    }

    stats.baseline_round_ns_per_query.reserve(static_cast<size_t>(timing_rounds));
    stats.tree_round_ns_per_query.reserve(static_cast<size_t>(timing_rounds));
    for (int i = 0; i < timing_rounds; ++i) {
        stats.baseline_round_ns_per_query.push_back(RunBaselineRound(scenario, grid, samples));
    }
    for (int i = 0; i < timing_rounds; ++i) {
        stats.tree_round_ns_per_query.push_back(RunTreeRound(scenario, tree, samples));
    }

    stats.baseline_mean_ns_per_query = Mean(stats.baseline_round_ns_per_query);
    stats.baseline_median_ns_per_query = Median(stats.baseline_round_ns_per_query);
    stats.tree_mean_ns_per_query = Mean(stats.tree_round_ns_per_query);
    stats.tree_median_ns_per_query = Median(stats.tree_round_ns_per_query);
    if (stats.tree_median_ns_per_query > 0.0) {
        stats.speedup_ratio = stats.baseline_median_ns_per_query / stats.tree_median_ns_per_query;
    }
    return stats;
}

ScenarioMetrics MakeScenarioMetrics(
    const ScenarioDefinition& scenario,
    const AccuracyStats& accuracy,
    const TimingStats& timing,
    uint64_t seed,
    size_t grid_steady_bytes,
    const AdaptiveQuadtreeSummary& summary,
    const AdaptiveQuadtreeMemoryFootprint& memory,
    double tree_build_ms
) {
    ScenarioMetrics metrics;
    metrics.scenario_name = scenario.name;
    metrics.sampling_mode = SamplingModeToString(scenario.sampling_mode);
    metrics.distribution = DistributionModeToString(scenario.distribution);
    metrics.sample_count = accuracy.sample_count;
    metrics.seed = seed;
    metrics.hit_count_baseline = accuracy.hit_count_baseline;
    metrics.hit_count_tree = accuracy.hit_count_tree;
    metrics.mae = accuracy.mae;
    metrics.rmse = accuracy.rmse;
    metrics.max_abs_error = accuracy.max_abs_error;
    metrics.p95_abs_error = accuracy.p95_abs_error;
    metrics.baseline_mean_ns_per_query = timing.baseline_mean_ns_per_query;
    metrics.baseline_median_ns_per_query = timing.baseline_median_ns_per_query;
    metrics.tree_mean_ns_per_query = timing.tree_mean_ns_per_query;
    metrics.tree_median_ns_per_query = timing.tree_median_ns_per_query;
    metrics.speedup_ratio = timing.speedup_ratio;
    metrics.grid_steady_bytes = grid_steady_bytes;
    metrics.tree_steady_bytes = memory.steady_bytes;
    metrics.tree_peak_incremental_bytes = memory.peak_build_incremental_bytes;
    metrics.tree_total_peak_bytes = grid_steady_bytes + memory.peak_build_incremental_bytes;
    metrics.tree_build_ms = tree_build_ms;
    metrics.leaf_count = summary.leaf_count;
    metrics.zero_leaf_count = summary.zero_leaf_count;
    metrics.smooth_leaf_count = summary.smooth_leaf_count;
    metrics.boundary_leaf_count = summary.boundary_leaf_count;
    return metrics;
}

ScenarioMetrics MakeSkippedScenarioMetrics(
    const ScenarioDefinition& scenario,
    uint64_t seed,
    size_t grid_steady_bytes,
    const AdaptiveQuadtreeSummary& summary,
    const AdaptiveQuadtreeMemoryFootprint& memory,
    double tree_build_ms,
    const std::string& reason
) {
    ScenarioMetrics metrics;
    metrics.scenario_name = scenario.name;
    metrics.sampling_mode = SamplingModeToString(scenario.sampling_mode);
    metrics.distribution = DistributionModeToString(scenario.distribution);
    metrics.seed = seed;
    metrics.grid_steady_bytes = grid_steady_bytes;
    metrics.tree_steady_bytes = memory.steady_bytes;
    metrics.tree_peak_incremental_bytes = memory.peak_build_incremental_bytes;
    metrics.tree_total_peak_bytes = grid_steady_bytes + memory.peak_build_incremental_bytes;
    metrics.tree_build_ms = tree_build_ms;
    metrics.leaf_count = summary.leaf_count;
    metrics.zero_leaf_count = summary.zero_leaf_count;
    metrics.smooth_leaf_count = summary.smooth_leaf_count;
    metrics.boundary_leaf_count = summary.boundary_leaf_count;
    metrics.skipped = true;
    metrics.skip_reason = reason;
    return metrics;
}

void WriteMetricsCsv(const std::string& output_path, const std::vector<ScenarioMetrics>& metrics) {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open benchmark metrics CSV: " + output_path);
    }

    file << "scenario_name,sampling_mode,distribution,sample_count,seed,hit_count_baseline,hit_count_tree,"
         << "mae,rmse,max_abs_error,p95_abs_error,"
         << "baseline_mean_ns_per_query,baseline_median_ns_per_query,"
         << "tree_mean_ns_per_query,tree_median_ns_per_query,speedup_ratio,"
         << "grid_steady_bytes,tree_steady_bytes,tree_peak_incremental_bytes,tree_total_peak_bytes,"
         << "tree_build_ms,leaf_count,zero_leaf_count,smooth_leaf_count,boundary_leaf_count,skipped,skip_reason\n";
    file << std::fixed << std::setprecision(6);
    for (const ScenarioMetrics& row : metrics) {
        file << row.scenario_name << ","
             << row.sampling_mode << ","
             << row.distribution << ","
             << row.sample_count << ","
             << row.seed << ","
             << row.hit_count_baseline << ","
             << row.hit_count_tree << ","
             << row.mae << ","
             << row.rmse << ","
             << row.max_abs_error << ","
             << row.p95_abs_error << ","
             << row.baseline_mean_ns_per_query << ","
             << row.baseline_median_ns_per_query << ","
             << row.tree_mean_ns_per_query << ","
             << row.tree_median_ns_per_query << ","
             << row.speedup_ratio << ","
             << row.grid_steady_bytes << ","
             << row.tree_steady_bytes << ","
             << row.tree_peak_incremental_bytes << ","
             << row.tree_total_peak_bytes << ","
             << row.tree_build_ms << ","
             << row.leaf_count << ","
             << row.zero_leaf_count << ","
             << row.smooth_leaf_count << ","
             << row.boundary_leaf_count << ","
             << (row.skipped ? 1 : 0) << ","
             << row.skip_reason << "\n";
    }
}

void WriteSummary(
    const std::string& output_path,
    const std::string& csv_path,
    const IndexRegion& region,
    const EnergyGrid2D& grid,
    const AdaptiveSplitConfig& split_config,
    const AdaptiveQuadtreeSummary& summary,
    const AdaptiveQuadtreeMemoryFootprint& memory,
    size_t grid_steady_bytes,
    uint64_t seed,
    int full_samples,
    int boundary_samples,
    int warmup_rounds,
    int timing_rounds,
    double tree_build_ms,
    const std::vector<ScenarioMetrics>& metrics
) {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open benchmark summary output: " + output_path);
    }

    file << std::fixed << std::setprecision(6);
    file << "input_csv=" << csv_path << "\n";
    file << "region_ix_min=" << region.ix_min << "\n";
    file << "region_ix_max=" << region.ix_max << "\n";
    file << "region_iy_min=" << region.iy_min << "\n";
    file << "region_iy_max=" << region.iy_max << "\n";
    file << "seed=" << seed << "\n";
    file << "full_samples=" << full_samples << "\n";
    file << "boundary_samples=" << boundary_samples << "\n";
    file << "warmup_rounds=" << warmup_rounds << "\n";
    file << "timing_rounds=" << timing_rounds << "\n";
    file << "grid_width=" << grid.width << "\n";
    file << "grid_height=" << grid.height << "\n";
    file << "grid_resolution_m=" << grid.resolution_m << "\n";
    file << "fixed_height_m=" << grid.fixed_height_m << "\n";
    file << "gradient_threshold_rel=" << split_config.gradient_threshold_rel << "\n";
    file << "center_error_threshold_abs=" << split_config.center_error_threshold_abs << "\n";
    file << "min_block_size=" << split_config.min_block_size << "\n";
    file << "zero_epsilon=" << split_config.zero_epsilon << "\n";
    file << "tree_build_ms=" << tree_build_ms << "\n";
    file << "leaf_count=" << summary.leaf_count << "\n";
    file << "zero_leaf_count=" << summary.zero_leaf_count << "\n";
    file << "smooth_leaf_count=" << summary.smooth_leaf_count << "\n";
    file << "boundary_leaf_count=" << summary.boundary_leaf_count << "\n";
    file << "grid_steady_bytes=" << grid_steady_bytes << "\n";
    file << "tree_steady_bytes=" << memory.steady_bytes << "\n";
    file << "tree_peak_incremental_bytes=" << memory.peak_build_incremental_bytes << "\n";
    file << "tree_total_peak_bytes=" << (grid_steady_bytes + memory.peak_build_incremental_bytes) << "\n";
    file << "tree_node_bytes=" << memory.node_bytes << "\n";
    file << "tree_leaf_bytes=" << memory.leaf_bytes << "\n";
    file << "tree_query_meta_bytes=" << memory.query_meta_bytes << "\n";
    file << "tree_prefix_peak_bytes=" << memory.prefix_peak_bytes << "\n";
    file << "baseline_build_ms=0.000000\n";

    for (const ScenarioMetrics& row : metrics) {
        file << "\n[" << row.scenario_name << "]\n";
        file << "sampling_mode=" << row.sampling_mode << "\n";
        file << "distribution=" << row.distribution << "\n";
        file << "sample_count=" << row.sample_count << "\n";
        file << "hit_count_baseline=" << row.hit_count_baseline << "\n";
        file << "hit_count_tree=" << row.hit_count_tree << "\n";
        file << "mae=" << row.mae << "\n";
        file << "rmse=" << row.rmse << "\n";
        file << "max_abs_error=" << row.max_abs_error << "\n";
        file << "p95_abs_error=" << row.p95_abs_error << "\n";
        file << "baseline_mean_ns_per_query=" << row.baseline_mean_ns_per_query << "\n";
        file << "baseline_median_ns_per_query=" << row.baseline_median_ns_per_query << "\n";
        file << "tree_mean_ns_per_query=" << row.tree_mean_ns_per_query << "\n";
        file << "tree_median_ns_per_query=" << row.tree_median_ns_per_query << "\n";
        file << "speedup_ratio=" << row.speedup_ratio << "\n";
        file << "skipped=" << (row.skipped ? 1 : 0) << "\n";
        file << "skip_reason=" << row.skip_reason << "\n";
    }
}

void PrintSummaryToStdout(
    const EnergyGrid2D& grid,
    const AdaptiveQuadtreeSummary& summary,
    const AdaptiveQuadtreeMemoryFootprint& memory,
    size_t grid_steady_bytes,
    double tree_build_ms,
    const std::vector<ScenarioMetrics>& metrics
) {
    std::cout << "Grid: " << grid.width << " x " << grid.height
              << ", resolution=" << grid.resolution_m << " m\n";
    std::cout << "Leaves: " << summary.leaf_count
              << " [zero=" << summary.zero_leaf_count
              << ", smooth=" << summary.smooth_leaf_count
              << ", boundary=" << summary.boundary_leaf_count << "]\n";
    std::cout << "Build: " << tree_build_ms << " ms\n";
    std::cout << "Memory: grid=" << grid_steady_bytes
              << " B, tree_steady=" << memory.steady_bytes
              << " B, tree_peak_incremental=" << memory.peak_build_incremental_bytes
              << " B\n";
    for (const ScenarioMetrics& row : metrics) {
        std::cout << row.scenario_name << ": ";
        if (row.skipped) {
            std::cout << "skipped (" << row.skip_reason << ")\n";
            continue;
        }
        std::cout << "MAE=" << row.mae
                  << ", RMSE=" << row.rmse
                  << ", baseline_median_ns/query=" << row.baseline_median_ns_per_query
                  << ", tree_median_ns/query=" << row.tree_median_ns_per_query
                  << ", speedup=" << row.speedup_ratio << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const ParsedArgs args = ParseArgs(argc, argv);
        if (argc == 1 || HasFlag(args, "--help")) {
            PrintUsage();
            return 0;
        }
        if (!HasOption(args, "--input-csv")) {
            throw std::runtime_error("--input-csv is required");
        }

        const std::string input_csv = GetString(args, "--input-csv", "");
        const std::optional<double> explicit_resolution = HasOption(args, "--grid-resolution")
            ? std::optional<double>(GetDouble(args, "--grid-resolution", 1.0))
            : std::nullopt;
        EnergyGrid2D grid = EnergyGrid2D::LoadFromCsv(input_csv, explicit_resolution);

        AdaptiveSplitConfig split_config;
        split_config.gradient_threshold_rel = GetDouble(args, "--gradient-threshold", split_config.gradient_threshold_rel);
        split_config.min_block_size = GetInt(args, "--min-block-size", split_config.min_block_size);
        if (HasOption(args, "--error-threshold")) {
            split_config.center_error_threshold_abs = GetDouble(
                args,
                "--error-threshold",
                split_config.center_error_threshold_abs
            );
        } else {
            split_config.center_error_threshold_abs = 0.05 * grid.max_positive_value(split_config.zero_epsilon);
        }

        const int full_samples = GetInt(args, "--full-samples", 100000);
        const int boundary_samples = GetInt(args, "--boundary-samples", 50000);
        const int warmup_rounds = GetInt(args, "--warmup-rounds", 3);
        const int timing_rounds = GetInt(args, "--timing-rounds", 20);
        const uint64_t seed = GetUInt64(args, "--seed", 20260419ULL);
        if (full_samples <= 0 || boundary_samples <= 0 || warmup_rounds < 0 || timing_rounds <= 0) {
            throw std::runtime_error("Sample counts and timing rounds must be positive; warmup rounds must be non-negative");
        }

        AdaptiveQuadtree2D tree;
        const auto build_start = std::chrono::steady_clock::now();
        tree.Build(grid, split_config);
        const auto build_end = std::chrono::steady_clock::now();
        const double tree_build_ms = std::chrono::duration<double, std::milli>(build_end - build_start).count();
        const size_t grid_steady_bytes = grid.EstimatedSteadyBytes();

        const std::vector<ScenarioDefinition> scenarios = {
            {"cell_uniform", SamplingMode::Cell, DistributionMode::Uniform, static_cast<size_t>(full_samples)},
            {"cell_boundary", SamplingMode::Cell, DistributionMode::Boundary, static_cast<size_t>(boundary_samples)},
            {"continuous_uniform", SamplingMode::Continuous, DistributionMode::Uniform, static_cast<size_t>(full_samples)},
            {"continuous_boundary", SamplingMode::Continuous, DistributionMode::Boundary, static_cast<size_t>(boundary_samples)},
        };

        const IndexRegion region = ResolveRegion(args, grid);
        std::vector<ScenarioMetrics> metrics;
        metrics.reserve(scenarios.size());
        if (!region.valid) {
            for (const ScenarioDefinition& scenario : scenarios) {
                metrics.push_back(MakeSkippedScenarioMetrics(
                    scenario,
                    seed,
                    grid_steady_bytes,
                    tree.summary(),
                    tree.memory_footprint(),
                    tree_build_ms,
                    region.invalid_reason
                ));
            }
        } else {
            const auto boundary_cells = BuildBoundaryCells(grid, region, split_config);
            std::mt19937_64 rng(seed);
            for (const ScenarioDefinition& scenario : scenarios) {
                const std::vector<QuerySample> samples = GenerateSamples(scenario, grid, region, boundary_cells, rng);
                if (scenario.distribution == DistributionMode::Boundary && samples.empty()) {
                    metrics.push_back(MakeSkippedScenarioMetrics(
                        scenario,
                        seed,
                        grid_steady_bytes,
                        tree.summary(),
                        tree.memory_footprint(),
                        tree_build_ms,
                        "empty_boundary_mask"
                    ));
                    continue;
                }

                const AccuracyStats accuracy = EvaluateAccuracy(scenario, grid, tree, samples);
                const TimingStats timing = EvaluateTiming(scenario, grid, tree, samples, warmup_rounds, timing_rounds);
                metrics.push_back(MakeScenarioMetrics(
                    scenario,
                    accuracy,
                    timing,
                    seed,
                    grid_steady_bytes,
                    tree.summary(),
                    tree.memory_footprint(),
                    tree_build_ms
                ));
            }
        }

        const std::string output_prefix = GetString(args, "--output-prefix", "fixed_height_quadtree_monte_carlo");
        WriteMetricsCsv(output_prefix + "_benchmark_metrics.csv", metrics);
        WriteSummary(
            output_prefix + "_benchmark_summary.txt",
            input_csv,
            region,
            grid,
            split_config,
            tree.summary(),
            tree.memory_footprint(),
            grid_steady_bytes,
            seed,
            full_samples,
            boundary_samples,
            warmup_rounds,
            timing_rounds,
            tree_build_ms,
            metrics
        );
        PrintSummaryToStdout(
            grid,
            tree.summary(),
            tree.memory_footprint(),
            grid_steady_bytes,
            tree_build_ms,
            metrics
        );
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "fixed_height_quadtree_monte_carlo error: " << ex.what() << std::endl;
        return 1;
    }
}
