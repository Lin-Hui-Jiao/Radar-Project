// vis_check.cpp — independent check of the visibility (0 / non-zero) of a vehicle_field512 field.
//
//   vis_check IDX.3idx PREFIX [N_RANDOM] [N_BOUNDARY] [SEED]
//
// field_export decides visibility with RTree::Intersect3d (its own segment / box traversal) and lineSegIntersectTri
// (Moller-Trumbore with EPSILON tolerances). This tool uses neither:
//   candidates : RTree::Search (plain box overlap) over 1 m pieces of the radar -> node segment;
//   test       : watertight segment / triangle test with signed volumes in long double. Closed triangles: a segment
//                through an edge or a vertex counts as a hit; a segment whose end point lies on the triangle plane
//                does not (the node itself may lie on a surface).
// Nodes: N_RANDOM uniform random nodes and N_BOUNDARY nodes on a shadow boundary (a face neighbour has the other
// state), where errors would show first. The column straight under the radar is skipped (zero by the kernel, not by
// occlusion). Geometry and radar are read from PREFIX.json (grid origin, rtree_offset, radar epsg2326 / abs_alt_m),
// with the same arithmetic as field_export. Exit code 0 = no disagreement.
#include "hiradar/occlusion_utils.h"
#include "json.hpp"
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <unordered_set>
#include <vector>

namespace {

typedef long double LD;

inline LD orient(const LD* a, const LD* b, const LD* c, const LD* d) {   // det[b - a, c - a, d - a]
    const LD bx = b[0] - a[0], by = b[1] - a[1], bz = b[2] - a[2];
    const LD cx = c[0] - a[0], cy = c[1] - a[1], cz = c[2] - a[2];
    const LD dx = d[0] - a[0], dy = d[1] - a[1], dz = d[2] - a[2];
    return bx * (cy * dz - cz * dy) - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx);
}

// 0 = miss, 1 = proper hit, 2 = hit through an edge or a vertex
int seg_tri(const LD* P, const LD* Q, const LD* A, const LD* B, const LD* C) {
    const LD sp = orient(A, B, C, P), sq = orient(A, B, C, Q);
    if (sp == 0 || sq == 0 || (sp > 0) == (sq > 0)) return 0;
    const LD o1 = orient(P, Q, A, B), o2 = orient(P, Q, B, C), o3 = orient(P, Q, C, A);
    if (!((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0))) return 0;
    return (o1 == 0 || o2 == 0 || o3 == 0) ? 2 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 6) { fprintf(stderr, "usage: vis_check IDX.3idx PREFIX [N_RANDOM] [N_BOUNDARY] [SEED]\n"); return 2; }
    const std::string pre = argv[2];
    const int nr = argc > 3 ? atoi(argv[3]) : 20000, nb = argc > 4 ? atoi(argv[4]) : 20000;
    const uint64_t seed = argc > 5 ? strtoull(argv[5], nullptr, 10) : 1;
    const pshc::json::Value js = pshc::json::parse_file(pre + ".json");
    const auto& g = js["grid"];
    const auto& pv = js["provenance"];
    const int nx = (int)g["nx"].as_int(), ny = (int)g["ny"].as_int(), nz = (int)g["nz"].as_int();
    const double E0 = g["origin"][0].as_num(), N0 = g["origin"][1].as_num(), Z0 = g["origin"][2].as_num();
    const double sx = g["dx"].as_num(), sy = g["dy"].as_num(), sz = g["dz"].as_num();
    if (g["crs"].as_str() != "EPSG:2326") { fprintf(stderr, "only metric (EPSG:2326) grids are supported\n"); return 2; }
    const double offx = pv["rtree_offset"][0].as_num(), offy = pv["rtree_offset"][1].as_num();
    const double RE = pv["radar"]["epsg2326"][0].as_num(), RN = pv["radar"]["epsg2326"][1].as_num();
    const double RZ = pv["radar"]["abs_alt_m"].as_num();

    std::unique_ptr<RTree3d> T(new RTree3d());
    if (!T->Load(argv[1])) { fprintf(stderr, "cannot load %s\n", argv[1]); return 1; }
    std::vector<float> F((size_t)nx * ny * nz);
    FILE* fp = fopen((pre + ".f32").c_str(), "rb");
    if (!fp || fread(F.data(), 4, F.size(), fp) != F.size()) { fprintf(stderr, "cannot read %s.f32\n", pre.c_str()); return 1; }
    fclose(fp);
    auto vis = [&](int i, int j, int k) { return F[((size_t)k * ny + j) * nx + i] > 0.f; };

    struct Node { int i, j, k, boundary; };
    std::vector<Node> nodes;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> Ui(0, nx - 1), Uj(0, ny - 1), Uk(0, nz - 1);
    std::uniform_int_distribution<int> Bi(1, nx - 2), Bj(1, ny - 2), Bk(1, nz - 2);
    while ((int)nodes.size() < nr) nodes.push_back({Ui(rng), Uj(rng), Uk(rng), 0});
    for (int got = 0; got < nb;) {
        const int i = Bi(rng), j = Bj(rng), k = Bk(rng);
        const bool v = vis(i, j, k);
        if (vis(i - 1, j, k) != v || vis(i + 1, j, k) != v || vis(i, j - 1, k) != v || vis(i, j + 1, k) != v ||
            vis(i, j, k - 1) != v || vis(i, j, k + 1) != v) { nodes.push_back({i, j, k, 1}); ++got; }
    }
    const double s[3] = {RE + offx, RN + offy, RZ};   // field_export: RX = r_E + off_x, RZ = ground + mount
    const LD P[3] = {s[0], s[1], s[2]};
    std::vector<int> res(nodes.size(), 0);   // -1 skipped, 0 agree, 1 field visible / test occluded, 2 the reverse
    size_t edge_hits = 0;
    #pragma omp parallel for schedule(dynamic, 64) reduction(+ : edge_hits)
    for (long q = 0; q < (long)nodes.size(); ++q) {
        const Node& nd = nodes[q];
        const double E = E0 + nd.i * sx, N = N0 + nd.j * sy;
        if (E == RE && N == RN) { res[q] = -1; continue; }
        const double e[3] = {E + offx, N + offy, Z0 + nd.k * sz};   // same arithmetic as field_export
        const LD Q[3] = {e[0], e[1], e[2]};
        const double len = std::sqrt((e[0] - s[0]) * (e[0] - s[0]) + (e[1] - s[1]) * (e[1] - s[1]) + (e[2] - s[2]) * (e[2] - s[2]));
        const int np = std::max(1, (int)std::ceil(len));
        std::unordered_set<const ValueType*> seen;
        int hit = 0;
        for (int p = 0; p < np && !hit; ++p) {
            double mn[3], mx[3];
            for (int a = 0; a < 3; ++a) {
                const double u = s[a] + (e[a] - s[a]) * ((double)p / np), w = s[a] + (e[a] - s[a]) * ((double)(p + 1) / np);
                mn[a] = std::min(u, w) - 1e-6; mx[a] = std::max(u, w) + 1e-6;
            }
            T->Search(mn, mx, [&](const ValueType& t) {
                if (!seen.insert(&t).second) return true;
                const LD A[3] = {std::get<0>(t), std::get<1>(t), std::get<2>(t)};
                const LD B[3] = {std::get<3>(t), std::get<4>(t), std::get<5>(t)};
                const LD C[3] = {std::get<6>(t), std::get<7>(t), std::get<8>(t)};
                hit = seg_tri(P, Q, A, B, C);
                return hit == 0;
            });
        }
        edge_hits += hit == 2;
        const bool fv = vis(nd.i, nd.j, nd.k);
        res[q] = fv == (hit == 0) ? 0 : fv ? 1 : 2;
    }
    size_t agree[2] = {0, 0}, d1[2] = {0, 0}, d2[2] = {0, 0}, skipped = 0, nvis[2] = {0, 0};
    for (size_t q = 0; q < nodes.size(); ++q) {
        const int b = nodes[q].boundary;
        if (res[q] < 0) { ++skipped; continue; }
        nvis[b] += vis(nodes[q].i, nodes[q].j, nodes[q].k);
        (res[q] == 0 ? agree : res[q] == 1 ? d1 : d2)[b]++;
    }
    printf("vis_check %s: radar E %.3f N %.3f alt %.6f m\n", pre.c_str(), RE, RN, RZ);
    for (int b = 0; b < 2; ++b)
        printf("  %s nodes (%zu visible in the field): agree %zu; field visible but occluded by the test %zu; field 0 but "
               "visible by the test %zu\n", b ? "boundary" : "random  ", nvis[b], agree[b], d1[b], d2[b]);
    printf("  skipped (column under the radar) %zu; hits found only through an edge or a vertex %zu\n", skipped, edge_hits);
    int shown = 0;
    for (size_t q = 0; q < nodes.size() && shown < 20; ++q)
        if (res[q] > 0) { printf("  disagree (i, j, k) = (%d, %d, %d)\n", nodes[q].i, nodes[q].j, nodes[q].k); ++shown; }
    return d1[0] + d1[1] + d2[0] + d2[1] == 0 ? 0 : 1;
}
