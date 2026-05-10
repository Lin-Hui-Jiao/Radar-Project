#include "hiradar/adaptive_quadtree_2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace hiradar {

namespace {

size_t PrefixIndex(int width, int x, int y) {
    return static_cast<size_t>(y) * static_cast<size_t>(width + 1) + static_cast<size_t>(x);
}

size_t VectorBytes(const std::vector<int>& values) {
    return values.capacity() * sizeof(int);
}

size_t VectorBytes(const std::vector<double>& values) {
    return values.capacity() * sizeof(double);
}

}  // namespace

AdaptiveQuadtree2D::AdaptiveQuadtree2D() = default;

const char* LeafStateToString(QuadtreeLeafState state) {
    switch (state) {
        case QuadtreeLeafState::Zero:
            return "zero";
        case QuadtreeLeafState::Smooth:
            return "smooth";
        case QuadtreeLeafState::Boundary:
            return "boundary";
        default:
            return "unknown";
    }
}

void AdaptiveQuadtree2D::Build(const EnergyGrid2D& grid, const AdaptiveSplitConfig& config) {
    if (grid.width <= 0 || grid.height <= 0 || grid.values.empty()) {
        throw std::runtime_error("AdaptiveQuadtree2D::Build requires a non-empty grid");
    }
    if (config.min_block_size <= 0) {
        throw std::runtime_error("min_block_size must be positive");
    }

    grid_ = &grid;
    root_index_ = -1;
    config_ = config;
    summary_ = AdaptiveQuadtreeSummary{};
    memory_footprint_ = AdaptiveQuadtreeMemoryFootprint{};
    query_metadata_ = QueryMetadata{};
    nodes_.clear();
    leaves_.clear();
    valid_prefix_.clear();
    positive_prefix_.clear();
    sum_prefix_.clear();

    query_metadata_.width = grid.width;
    query_metadata_.height = grid.height;
    query_metadata_.resolution_m = grid.resolution_m;
    query_metadata_.fixed_height_m = grid.fixed_height_m;
    query_metadata_.x_min_m = grid.x_min_m;
    query_metadata_.y_min_m = grid.y_min_m;
    query_metadata_.origin_frame = grid.origin_frame;

    summary_.original_width = grid.width;
    summary_.original_height = grid.height;
    summary_.original_points = grid.point_count();
    summary_.valid_points = grid.point_count();
    summary_.global_max_positive = grid.max_positive_value(config.zero_epsilon);
    if (config_.center_error_threshold_abs < 0.0) {
        config_.center_error_threshold_abs = 0.05 * summary_.global_max_positive;
    }

    summary_.padded_size = NextPowerOfTwo(std::max(grid.width, grid.height));
    const int prefix_width = grid.width + 1;
    const int prefix_height = grid.height + 1;
    valid_prefix_.assign(static_cast<size_t>(prefix_width) * prefix_height, 0);
    positive_prefix_.assign(static_cast<size_t>(prefix_width) * prefix_height, 0);
    sum_prefix_.assign(static_cast<size_t>(prefix_width) * prefix_height, 0.0);
    UpdatePeakMemoryEstimate();

    for (int y = 0; y < grid.height; ++y) {
        for (int x = 0; x < grid.width; ++x) {
            const size_t src_index = static_cast<size_t>(y) * grid.width + x;
            const int valid = 1;
            const int positive = grid.values[src_index] > config.zero_epsilon ? 1 : 0;
            const double sum = static_cast<double>(grid.values[src_index]);
            const size_t dst = PrefixIndex(grid.width, x + 1, y + 1);
            const size_t left = PrefixIndex(grid.width, x, y + 1);
            const size_t up = PrefixIndex(grid.width, x + 1, y);
            const size_t diag = PrefixIndex(grid.width, x, y);

            valid_prefix_[dst] = valid + valid_prefix_[left] + valid_prefix_[up] - valid_prefix_[diag];
            positive_prefix_[dst] = positive + positive_prefix_[left] + positive_prefix_[up] - positive_prefix_[diag];
            sum_prefix_[dst] = sum + sum_prefix_[left] + sum_prefix_[up] - sum_prefix_[diag];
        }
    }

    root_index_ = BuildNode(0, 0, summary_.padded_size, false);
    if (root_index_ < 0) {
        throw std::runtime_error("AdaptiveQuadtree2D failed to build a valid root node");
    }

    summary_.leaf_count = leaves_.size();

    valid_prefix_.clear();
    positive_prefix_.clear();
    sum_prefix_.clear();
    valid_prefix_.shrink_to_fit();
    positive_prefix_.shrink_to_fit();
    sum_prefix_.shrink_to_fit();
    grid_ = nullptr;
    UpdatePeakMemoryEstimate();
}

std::optional<float> AdaptiveQuadtree2D::QueryByIndex(int ix, int iy) const {
    if (root_index_ < 0) {
        return std::nullopt;
    }
    if (ix < 0 || ix >= query_metadata_.width || iy < 0 || iy >= query_metadata_.height) {
        return std::nullopt;
    }
    return QueryNode(root_index_, ix, iy);
}

std::optional<float> AdaptiveQuadtree2D::QueryByLocalMeters(double x_m, double y_m) const {
    if (root_index_ < 0 || query_metadata_.resolution_m <= 0.0) {
        return std::nullopt;
    }
    const int ix = static_cast<int>(std::floor((x_m - query_metadata_.x_min_m) / query_metadata_.resolution_m));
    const int iy = static_cast<int>(std::floor((y_m - query_metadata_.y_min_m) / query_metadata_.resolution_m));
    return QueryByIndex(ix, iy);
}

std::optional<float> AdaptiveQuadtree2D::QueryByGeodetic(double lon_deg, double lat_deg) const {
    if (root_index_ < 0 || !query_metadata_.origin_frame.IsInitialized()) {
        return std::nullopt;
    }
    const Vec3d enu = query_metadata_.origin_frame.ToENU(lon_deg, lat_deg, query_metadata_.fixed_height_m);
    return QueryByLocalMeters(enu.x, enu.y);
}

void AdaptiveQuadtree2D::ExportLeavesCsv(const std::string& output_path) const {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open leaf CSV output: " + output_path);
    }

    file << "x0,y0,size,value,state\n";
    for (const QuadtreeLeaf& leaf : leaves_) {
        file << leaf.x0 << ","
             << leaf.y0 << ","
             << leaf.size << ","
             << leaf.value << ","
             << LeafStateToString(static_cast<QuadtreeLeafState>(leaf.state)) << "\n";
    }
}

int AdaptiveQuadtree2D::BuildNode(int x0, int y0, int size, bool boundary_context) {
    const BlockStats stats = EvaluateBlock(x0, y0, size);
    if (stats.valid_count == 0) {
        return -1;
    }

    const bool is_zero_block = stats.positive_count == 0;
    const bool block_mixed = stats.positive_count > 0 && stats.positive_count < stats.valid_count;
    const bool reached_min_block = size <= config_.min_block_size;
    const bool should_split =
        !reached_min_block &&
        !is_zero_block &&
        (block_mixed ||
         stats.zero_nonzero_mixed ||
         stats.gradient_rel > config_.gradient_threshold_rel ||
         stats.center_error > config_.center_error_threshold_abs);

    Node node;
    node.x0 = x0;
    node.y0 = y0;
    node.size = size;

    const bool child_boundary_context = boundary_context || block_mixed || stats.zero_nonzero_mixed;
    if (should_split) {
        const int half = size / 2;
        const int child_coords[4][2] = {
            {x0, y0},
            {x0 + half, y0},
            {x0, y0 + half},
            {x0 + half, y0 + half},
        };

        int valid_child_count = 0;
        const int node_index = static_cast<int>(nodes_.size());
        nodes_.push_back(node);
        UpdatePeakMemoryEstimate();
        for (int child = 0; child < 4; ++child) {
            const int child_index = BuildNode(
                child_coords[child][0],
                child_coords[child][1],
                half,
                child_boundary_context
            );
            nodes_[node_index].children[child] = child_index;
            if (child_index >= 0) {
                ++valid_child_count;
            }
        }

        if (valid_child_count > 0) {
            return node_index;
        }
        nodes_.pop_back();
        UpdatePeakMemoryEstimate();
    }

    node.is_leaf = true;
    QuadtreeLeaf leaf;
    leaf.x0 = x0;
    leaf.y0 = y0;
    leaf.size = size;
    leaf.value = is_zero_block ? 0.0f : stats.mean;
    if (is_zero_block) {
        leaf.state = static_cast<uint8_t>(QuadtreeLeafState::Zero);
        ++summary_.zero_leaf_count;
    } else if (child_boundary_context) {
        leaf.state = static_cast<uint8_t>(QuadtreeLeafState::Boundary);
        ++summary_.boundary_leaf_count;
    } else {
        leaf.state = static_cast<uint8_t>(QuadtreeLeafState::Smooth);
        ++summary_.smooth_leaf_count;
    }

    node.leaf_index = static_cast<int>(leaves_.size());
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back(node);
    leaves_.push_back(leaf);
    UpdatePeakMemoryEstimate();
    return node_index;
}

AdaptiveQuadtree2D::BlockStats AdaptiveQuadtree2D::EvaluateBlock(int x0, int y0, int size) const {
    BlockStats stats;
    stats.x0 = x0;
    stats.y0 = y0;
    stats.size = size;

    const int x1 = std::max(0, std::min(x0, grid_->width));
    const int y1 = std::max(0, std::min(y0, grid_->height));
    const int x2 = std::max(0, std::min(x0 + size, grid_->width));
    const int y2 = std::max(0, std::min(y0 + size, grid_->height));

    stats.valid_count = RangeSumInt(valid_prefix_, x1, y1, x2, y2);
    stats.positive_count = RangeSumInt(positive_prefix_, x1, y1, x2, y2);
    stats.sum = RangeSumDouble(sum_prefix_, x1, y1, x2, y2);
    if (stats.valid_count == 0) {
        return stats;
    }
    stats.mean = static_cast<float>(stats.sum / static_cast<double>(stats.valid_count));

    const int corner_x_min = std::min(x0, grid_->width - 1);
    const int corner_y_min = std::min(y0, grid_->height - 1);
    const int corner_x_max = std::min(x0 + size - 1, grid_->width - 1);
    const int corner_y_max = std::min(y0 + size - 1, grid_->height - 1);
    const int center_x = std::min(x0 + size / 2, grid_->width - 1);
    const int center_y = std::min(y0 + size / 2, grid_->height - 1);

    stats.corners[0] = grid_->get(corner_x_min, corner_y_min);
    stats.corners[1] = grid_->get(corner_x_max, corner_y_min);
    stats.corners[2] = grid_->get(corner_x_min, corner_y_max);
    stats.corners[3] = grid_->get(corner_x_max, corner_y_max);
    stats.center = grid_->get(center_x, center_y);

    float max_corner = stats.corners[0];
    float min_corner = stats.corners[0];
    int zero_count = 0;
    for (float corner : stats.corners) {
        max_corner = std::max(max_corner, corner);
        min_corner = std::min(min_corner, corner);
        if (corner <= config_.zero_epsilon) {
            ++zero_count;
        }
    }

    stats.zero_nonzero_mixed = zero_count > 0 && zero_count < 4;
    const double denominator = std::max(static_cast<double>(max_corner), config_.zero_epsilon);
    stats.gradient_rel = static_cast<float>((max_corner - min_corner) / denominator);
    const float corner_mean = 0.25f * (stats.corners[0] + stats.corners[1] + stats.corners[2] + stats.corners[3]);
    stats.center_error = std::abs(stats.center - corner_mean);
    return stats;
}

std::optional<float> AdaptiveQuadtree2D::QueryNode(int node_index, int ix, int iy) const {
    if (node_index < 0 || node_index >= static_cast<int>(nodes_.size())) {
        return std::nullopt;
    }

    const Node& node = nodes_[node_index];
    if (node.is_leaf) {
        return leaves_[node.leaf_index].value;
    }

    const int half = node.size / 2;
    const int child_x = ix >= node.x0 + half ? 1 : 0;
    const int child_y = iy >= node.y0 + half ? 1 : 0;
    const int child_slot = child_y * 2 + child_x;
    const int child_index = node.children[child_slot];
    if (child_index < 0) {
        return std::nullopt;
    }
    return QueryNode(child_index, ix, iy);
}

void AdaptiveQuadtree2D::UpdatePeakMemoryEstimate() {
    memory_footprint_.node_bytes = nodes_.capacity() * sizeof(Node);
    memory_footprint_.leaf_bytes = leaves_.capacity() * sizeof(QuadtreeLeaf);
    memory_footprint_.query_meta_bytes = sizeof(QueryMetadata);

    const size_t current_prefix_bytes =
        VectorBytes(valid_prefix_) +
        VectorBytes(positive_prefix_) +
        VectorBytes(sum_prefix_);

    memory_footprint_.prefix_peak_bytes = std::max(memory_footprint_.prefix_peak_bytes, current_prefix_bytes);
    memory_footprint_.steady_bytes = sizeof(AdaptiveQuadtree2D) + memory_footprint_.node_bytes + memory_footprint_.leaf_bytes;
    memory_footprint_.peak_build_incremental_bytes = std::max(
        memory_footprint_.peak_build_incremental_bytes,
        memory_footprint_.steady_bytes + current_prefix_bytes
    );
}

int AdaptiveQuadtree2D::RangeSumInt(const std::vector<int>& prefix, int x0, int y0, int x1, int y1) const {
    if (x0 >= x1 || y0 >= y1) {
        return 0;
    }
    const size_t br = PrefixIndex(grid_->width, x1, y1);
    const size_t bl = PrefixIndex(grid_->width, x0, y1);
    const size_t tr = PrefixIndex(grid_->width, x1, y0);
    const size_t tl = PrefixIndex(grid_->width, x0, y0);
    return prefix[br] - prefix[bl] - prefix[tr] + prefix[tl];
}

double AdaptiveQuadtree2D::RangeSumDouble(const std::vector<double>& prefix, int x0, int y0, int x1, int y1) const {
    if (x0 >= x1 || y0 >= y1) {
        return 0.0;
    }
    const size_t br = PrefixIndex(grid_->width, x1, y1);
    const size_t bl = PrefixIndex(grid_->width, x0, y1);
    const size_t tr = PrefixIndex(grid_->width, x1, y0);
    const size_t tl = PrefixIndex(grid_->width, x0, y0);
    return prefix[br] - prefix[bl] - prefix[tr] + prefix[tl];
}

int AdaptiveQuadtree2D::NextPowerOfTwo(int value) const {
    int power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
}

}  // namespace hiradar
