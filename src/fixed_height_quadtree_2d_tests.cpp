#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "hiradar/adaptive_quadtree_2d.hpp"

namespace {

using hiradar::AdaptiveQuadtree2D;
using hiradar::AdaptiveSplitConfig;
using hiradar::EnergyGrid2D;
using hiradar::LocalTangentFrame;
using hiradar::QuadtreeLeafState;

void Expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

EnergyGrid2D MakeGrid(int width, int height, const std::vector<float>& values) {
    LocalTangentFrame frame(114.1670, 22.2806, 50.0);
    return EnergyGrid2D::CreateForTesting(width, height, 1.0, 50.0, frame, 0.0, 0.0, values);
}

void TestRoundTrip() {
    LocalTangentFrame frame(114.1670, 22.2806, 50.0);
    const Vec3d enu = frame.ToENU(114.1675, 22.2811, 63.5);
    const Vec3d geo = frame.FromENU(enu.x, enu.y, enu.z);
    const Vec3d roundtrip_enu = frame.ToENU(geo.x, geo.y, geo.z);
    const double error = std::sqrt(
        (roundtrip_enu.x - enu.x) * (roundtrip_enu.x - enu.x) +
        (roundtrip_enu.y - enu.y) * (roundtrip_enu.y - enu.y) +
        (roundtrip_enu.z - enu.z) * (roundtrip_enu.z - enu.z)
    );
    Expect(error < 0.01, "LocalTangentFrame round-trip exceeds 1 cm");
}

void TestZeroField() {
    EnergyGrid2D grid = MakeGrid(8, 8, std::vector<float>(64, 0.0f));
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});
    Expect(tree.leaves().size() == 1, "Zero field should collapse to one leaf");
    Expect(tree.leaves().front().state == static_cast<uint8_t>(QuadtreeLeafState::Zero), "Zero field leaf state mismatch");
}

void TestConstantField() {
    EnergyGrid2D grid = MakeGrid(8, 8, std::vector<float>(64, 3.0f));
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});
    Expect(tree.leaves().size() == 1, "Constant field should collapse to one leaf");
    Expect(tree.leaves().front().state == static_cast<uint8_t>(QuadtreeLeafState::Smooth), "Constant field leaf state mismatch");
}

void TestBoundaryField() {
    std::vector<float> values(64, 0.0f);
    for (int y = 0; y < 8; ++y) {
        for (int x = 4; x < 8; ++x) {
            values[static_cast<size_t>(y) * 8 + x] = 2.0f;
        }
    }
    EnergyGrid2D grid = MakeGrid(8, 8, values);
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});
    bool has_boundary = false;
    for (const auto& leaf : tree.leaves()) {
        if (leaf.state == static_cast<uint8_t>(QuadtreeLeafState::Boundary)) {
            has_boundary = true;
            break;
        }
    }
    Expect(has_boundary, "Boundary field should produce boundary leaves");
}

void TestRectangularPadding() {
    EnergyGrid2D grid = MakeGrid(500, 503, std::vector<float>(500 * 503, 0.0f));
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});
    Expect(tree.summary().padded_size == 512, "Rectangular padding should expand to 512");
    Expect(tree.leaves().size() == 1, "Zero rectangular field should not create padding artifacts");
}

void TestQueryConsistency() {
    std::vector<float> values(16, 0.0f);
    values[5] = 1.5f;
    values[6] = 1.5f;
    values[9] = 1.5f;
    values[10] = 1.5f;
    EnergyGrid2D grid = MakeGrid(4, 4, values);
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});

    const auto by_index = tree.QueryByIndex(1, 1);
    const auto by_local = tree.QueryByLocalMeters(1.2, 1.2);
    const Vec3d geo = grid.origin_frame.FromENU(1.2, 1.2, 0.0);
    const auto by_geo = tree.QueryByGeodetic(geo.x, geo.y);

    Expect(by_index.has_value(), "Index query should hit");
    Expect(by_local.has_value(), "Local-meter query should hit");
    Expect(by_geo.has_value(), "Geodetic query should hit");
    Expect(std::abs(*by_index - *by_local) < 1e-6f, "Index/local query mismatch");
    Expect(std::abs(*by_index - *by_geo) < 1e-6f, "Index/geodetic query mismatch");
    Expect(!tree.QueryByIndex(10, 10).has_value(), "Out-of-bounds query should miss");
}

void TestInteriorHotspotForcesSplit() {
    std::vector<float> values(16, 0.0f);
    values[5] = 1.0f;
    EnergyGrid2D grid = MakeGrid(4, 4, values);
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});

    Expect(tree.leaves().size() > 1, "Interior hotspot should not collapse into one smooth leaf");
    const auto hotspot = tree.QueryByIndex(1, 1);
    const auto background = tree.QueryByIndex(3, 3);
    Expect(hotspot.has_value(), "Hotspot query should hit");
    Expect(background.has_value(), "Background query should hit");
    Expect(std::abs(*hotspot - 1.0f) < 1e-6f, "Hotspot should be isolated down to its cell value");
    Expect(std::abs(*background) < 1e-6f, "Background cell should remain zero");
}

void TestDetachedQueriesAndMemoryFootprint() {
    std::vector<float> values(16, 0.0f);
    values[5] = 2.0f;
    values[6] = 2.0f;
    values[9] = 2.0f;
    values[10] = 2.0f;
    EnergyGrid2D grid = MakeGrid(4, 4, values);
    AdaptiveQuadtree2D tree;
    tree.Build(grid, AdaptiveSplitConfig{});

    const size_t grid_bytes = grid.EstimatedSteadyBytes();
    const auto memory = tree.memory_footprint();
    Expect(grid_bytes > 0, "Grid steady-byte estimate should be positive");
    Expect(memory.steady_bytes > 0, "Tree steady-byte estimate should be positive");
    Expect(memory.peak_build_incremental_bytes >= memory.steady_bytes, "Tree peak build bytes should cover steady bytes");
    Expect(memory.query_meta_bytes > 0, "Tree query metadata bytes should be positive");

    grid.values.clear();
    grid.values.shrink_to_fit();
    const auto value = tree.QueryByIndex(1, 1);
    Expect(value.has_value(), "Detached tree query should still hit after source grid values are released");
    Expect(std::abs(*value - 2.0f) < 1e-6f, "Detached tree query returned an unexpected value");
}

void TestOptionalCsvLoad(const std::optional<std::string>& csv_path) {
    if (!csv_path.has_value()) {
        return;
    }
    EnergyGrid2D grid = EnergyGrid2D::LoadFromCsv(*csv_path);
    Expect(grid.width > 0 && grid.height > 0, "CSV grid dimensions should be positive");
    AdaptiveQuadtree2D tree;
    AdaptiveSplitConfig cfg;
    cfg.center_error_threshold_abs = 0.05 * grid.max_positive_value(cfg.zero_epsilon);
    tree.Build(grid, cfg);
    Expect(tree.summary().leaf_count > 0, "CSV-backed quadtree should have leaves");
    Expect(tree.memory_footprint().steady_bytes > 0, "CSV-backed quadtree should report steady memory");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        TestRoundTrip();
        TestZeroField();
        TestConstantField();
        TestBoundaryField();
        TestRectangularPadding();
        TestQueryConsistency();
        TestInteriorHotspotForcesSplit();
        TestDetachedQueriesAndMemoryFootprint();
        if (argc > 1) {
            TestOptionalCsvLoad(std::string(argv[1]));
        }
        std::cout << "fixed_height_quadtree_2d_tests: all checks passed" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "fixed_height_quadtree_2d_tests failed: " << ex.what() << std::endl;
        return 1;
    }
}
