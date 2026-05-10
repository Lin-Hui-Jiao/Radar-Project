#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "hiradar/energy_grid_2d.hpp"

namespace hiradar {

enum class QuadtreeLeafState : uint8_t {
    Zero = 0,
    Smooth = 1,
    Boundary = 2,
};

struct AdaptiveSplitConfig {
    int min_block_size = 1;
    double gradient_threshold_rel = 0.15;
    double center_error_threshold_abs = -1.0;
    double zero_epsilon = 1e-12;
};

struct QuadtreeLeaf {
    int x0 = 0;
    int y0 = 0;
    int size = 0;
    float value = 0.0f;
    uint8_t state = static_cast<uint8_t>(QuadtreeLeafState::Zero);
};

struct AdaptiveQuadtreeSummary {
    int original_width = 0;
    int original_height = 0;
    int padded_size = 0;
    size_t original_points = 0;
    size_t valid_points = 0;
    size_t leaf_count = 0;
    size_t zero_leaf_count = 0;
    size_t smooth_leaf_count = 0;
    size_t boundary_leaf_count = 0;
    double global_max_positive = 0.0;
};

struct AdaptiveQuadtreeMemoryFootprint {
    size_t steady_bytes = 0;
    size_t peak_build_incremental_bytes = 0;
    size_t node_bytes = 0;
    size_t leaf_bytes = 0;
    size_t query_meta_bytes = 0;
    size_t prefix_peak_bytes = 0;
};

class AdaptiveQuadtree2D {
public:
    AdaptiveQuadtree2D();

    void Build(const EnergyGrid2D& grid, const AdaptiveSplitConfig& config);

    const std::vector<QuadtreeLeaf>& leaves() const { return leaves_; }
    const AdaptiveQuadtreeSummary& summary() const { return summary_; }
    const AdaptiveQuadtreeMemoryFootprint& memory_footprint() const { return memory_footprint_; }

    std::optional<float> QueryByIndex(int ix, int iy) const;
    std::optional<float> QueryByLocalMeters(double x_m, double y_m) const;
    std::optional<float> QueryByGeodetic(double lon_deg, double lat_deg) const;

    void ExportLeavesCsv(const std::string& output_path) const;

private:
    struct Node {
        int x0 = 0;
        int y0 = 0;
        int size = 0;
        bool is_leaf = false;
        int leaf_index = -1;
        int children[4] = {-1, -1, -1, -1};
    };

    struct BlockStats {
        int x0 = 0;
        int y0 = 0;
        int size = 0;
        int valid_count = 0;
        int positive_count = 0;
        double sum = 0.0;
        float corners[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float center = 0.0f;
        float mean = 0.0f;
        float gradient_rel = 0.0f;
        float center_error = 0.0f;
        bool zero_nonzero_mixed = false;
    };

    struct QueryMetadata {
        int width = 0;
        int height = 0;
        double resolution_m = 1.0;
        double fixed_height_m = 0.0;
        double x_min_m = 0.0;
        double y_min_m = 0.0;
        LocalTangentFrame origin_frame;
    };

    int BuildNode(int x0, int y0, int size, bool boundary_context);
    BlockStats EvaluateBlock(int x0, int y0, int size) const;
    std::optional<float> QueryNode(int node_index, int ix, int iy) const;
    void UpdatePeakMemoryEstimate();

    int RangeSumInt(const std::vector<int>& prefix, int x0, int y0, int x1, int y1) const;
    double RangeSumDouble(const std::vector<double>& prefix, int x0, int y0, int x1, int y1) const;
    int NextPowerOfTwo(int value) const;

    const EnergyGrid2D* grid_ = nullptr;
    int root_index_ = -1;
    AdaptiveSplitConfig config_;
    AdaptiveQuadtreeSummary summary_;
    AdaptiveQuadtreeMemoryFootprint memory_footprint_;
    QueryMetadata query_metadata_;
    std::vector<Node> nodes_;
    std::vector<QuadtreeLeaf> leaves_;
    std::vector<int> valid_prefix_;
    std::vector<int> positive_prefix_;
    std::vector<double> sum_prefix_;
};

const char* LeafStateToString(QuadtreeLeafState state);

}  // namespace hiradar
