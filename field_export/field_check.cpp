// field_check.cpp — acceptance checks for field_export output (<prefix>.f32 + <prefix>.json).
//   field_check stats  PREFIX                       sha256, visible fraction, column transition histogram
//   field_check csv    PREFIX CSV [TOL]             compare one layer with a fixed_height_result_*.csv where both are
//                                                   non-zero (the CSV uses DEM occlusion; amplitudes must agree)
//   field_check random PREFIX INDEX [N] [SEED]      recompute N random nodes with Intersect3d + the amplitude kernel,
//                                                   coordinates derived here independently of field_export
#include "hiradar/occlusion_utils.h"
#include "radar_kernels.hpp"
#include "sha256.hpp"
#include "json.hpp"
#include <chrono>
#include <fstream>
#include <memory>
#include <random>

bool IntersectCallback(ValueType id, const double a[3], const double b[3]);
double HeightCallback(ValueType id, const double x, const double y);

namespace {

struct Field {
    pshc::json::Value js;
    int nx = 0, ny = 0, nz = 0;
    double dx = 0, dy = 0, dz = 0, o[3] = {0, 0, 0};
    bool lonlat = false;
    std::vector<float> v;
    float at(int i, int j, int k) const { return v[((size_t)k * ny + j) * nx + i]; }
};

Field load(const std::string& prefix, bool data = true) {
    Field f;
    f.js = pshc::json::parse_file(prefix + ".json");
    const auto& g = f.js["grid"];
    f.nx = (int)g["nx"].as_int(); f.ny = (int)g["ny"].as_int(); f.nz = (int)g["nz"].as_int();
    f.dx = g["dx"].as_num(); f.dy = g["dy"].as_num(); f.dz = g["dz"].as_num();
    for (int a = 0; a < 3; ++a) f.o[a] = g["origin"][a].as_num();
    f.lonlat = g["crs"].as_str() == "EPSG:4326";
    if (data) {
        size_t n = (size_t)f.nx * f.ny * f.nz; f.v.resize(n);
        FILE* fp = fopen((prefix + ".f32").c_str(), "rb");
        if (!fp || fread(f.v.data(), 4, n, fp) != n) { fprintf(stderr, "cannot read %s.f32\n", prefix.c_str()); exit(1); }
        fclose(fp);
    }
    return f;
}

int stats(const std::string& prefix) {
    Field f = load(prefix);
    pshc::Sha256 sh; sh.update(f.v.data(), f.v.size() * 4);
    const std::string d = sh.hex(), want = f.js["sha256"].as_str();
    size_t nvis = 0, hist[8] = {0}, nonmono = 0;
    const size_t plane = (size_t)f.nx * f.ny;
    for (size_t c = 0; c < plane; ++c) {
        int tr = 0; bool down = false; bool prev = false;
        for (int k = 0; k < f.nz; ++k) {
            bool b = f.v[k * plane + c] > 0.f; nvis += b;
            if (k && b != prev) { ++tr; if (prev && !b) down = true; }
            prev = b;
        }
        hist[std::min(tr, 7)]++; nonmono += down;
    }
    printf("%s: %d x %d x %d, sha256 %s (%s)\n", prefix.c_str(), f.nx, f.ny, f.nz, d.c_str(), d == want ? "matches json" : "MISMATCH");
    printf("visible (non-zero) fraction %.6f (%zu of %zu)\n", (double)nvis / f.v.size(), nvis, f.v.size());
    printf("columns by number of transitions along z:");
    for (int t = 0; t < 8; ++t) printf(" %d%s:%zu", t, t == 7 ? "+" : "", hist[t]);
    printf("\ncolumns with non-zero -> zero going up (true 3D occlusion): %zu of %zu\n", nonmono, plane);
    return d == want ? 0 : 1;
}

int csv(const std::string& prefix, const std::string& path, double tol) {
    Field f = load(prefix);
    std::ifstream in(path);
    if (!in) { fprintf(stderr, "cannot open %s\n", path.c_str()); return 1; }
    std::string line; std::getline(in, line);
    size_t row = 0, both = 0, bad = 0, axis_bad = 0, only_csv = 0, only_field = 0, kbad = 0;
    double maxd = 0, maxrel = 0; int klayer = -1;
    while (std::getline(in, line)) {
        double c[7]; if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6]) != 7) continue;
        const int i = (int)(row % f.nx), j = (int)(row / f.nx); ++row;
        if (j >= f.ny) { fprintf(stderr, "CSV has more rows than the grid\n"); return 1; }
        const int k = (int)std::lround((c[2] - f.o[2]) / f.dz);
        if (klayer < 0) klayer = k;
        if (k != klayer || k < 0 || k >= f.nz) { ++kbad; continue; }
        if (std::fabs(f.o[0] + i * f.dx - c[0]) > 6e-11 || std::fabs(f.o[1] + j * f.dy - c[1]) > 6e-11) ++axis_bad;
        const double a = c[6], b = f.at(i, j, k);
        if (a != 0 && b != 0) {
            ++both; double dd = std::fabs(a - b); maxd = std::max(maxd, dd); maxrel = std::max(maxrel, dd / a);
            if (dd > tol) ++bad;
        } else if (a != 0) ++only_csv; else if (b != 0) ++only_field;
    }
    printf("%s vs %s: layer k=%d (alt %.3f m), %zu rows; lon/lat axis mismatches %zu; rows off-layer %zu\n", prefix.c_str(),
           path.c_str(), klayer, f.o[2] + klayer * f.dz, row, axis_bad, kbad);
    printf("both non-zero %zu: max |diff| %.3g W/m^2 (max rel %.3g), above tolerance %.3g: %zu; non-zero only in CSV %zu, only in field %zu\n",
           both, maxd, maxrel, tol, bad, only_csv, only_field);
    return bad == 0 && axis_bad == 0 && row == (size_t)f.nx * f.ny ? 0 : 1;
}

int random_check(const std::string& prefix, const std::string& index, size_t n, uint64_t seed) {
    Field f = load(prefix);
    const auto& pv = f.js["provenance"];
    const auto& rd = pv["radar"];
    const bool vfloat = pv["data_version"].as_str() == "V-float";
    const double offx = pv["rtree_offset"][0].as_num(), offy = pv["rtree_offset"][1].as_num();
    auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<RTree3d> T(new RTree3d());
    if (!T->Load(index.c_str())) { fprintf(stderr, "cannot load %s\n", index.c_str()); return 1; }
    printf("index loaded in %.1f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    PJ_CONTEXT* ctx = proj_context_create();
    PJ* pj = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:2326", NULL);

    // radar: re-derive the occlusion position and the ground height from the inputs
    const bool in_ll = rd["input_crs"].as_str() == "EPSG:4326";
    const double ix = rd["input"][0].as_num(), iy = rd["input"][1].as_num();
    double rlon, rlat, rE, rN;
    if (in_ll) {
        rlon = ix; rlat = iy; double olon = ix, olat = iy;
        if (vfloat) {
            float a[3], b[3]; ReplayKernel::DecimalToDMS(a, (float)ix); ReplayKernel::DecimalToDMS(b, (float)iy);
            olon = DMSToDecimal(a); olat = DMSToDecimal(b);          // src/occlusion_utils.cpp
        }
        PJ_COORD c = proj_trans(pj, PJ_FWD, proj_coord(olat, olon, 0, 0)); rE = c.xy.y; rN = c.xy.x;
    } else {
        rE = ix; rN = iy;
        PJ_COORD c = proj_trans(pj, PJ_INV, proj_coord(rN, rE, 0, 0)); rlat = c.lp.lam; rlon = c.lp.phi;
    }
    const double RX = rE + offx, RY = rN + offy;
    const double ground = T->Getheight3d(RX, RY, HeightCallback), RZ = ground + rd["mount_height_m"].as_num();
    printf("radar: local (%.4f, %.4f) ground %.6f abs alt %.6f m; json says local (%.4f, %.4f) abs alt %.6f m\n", RX, RY,
           ground, RZ, rd["rtree_local"][0].as_num(), rd["rtree_local"][1].as_num(), rd["abs_alt_m"].as_num());
    Meta km; km.rlon = rlon; km.rlat = rlat; km.ralt = vfloat ? (double)(float)RZ : RZ;
    km.minphi = rd["elevation_range_deg"][0].as_num(); km.maxphi = rd["elevation_range_deg"][1].as_num();
    ReplayKernel KR; CleanKernel KC;
    if (vfloat) KR.init(km); else KC.init(km);

    std::mt19937_64 rng(seed);
    size_t mism = 0, nvis = 0;
    for (size_t t = 0; t < n; ++t) {
        const int i = (int)(rng() % f.nx), j = (int)(rng() % f.ny), k = (int)(rng() % f.nz);
        double lon, lat, E, N;
        if (f.lonlat) {
            lon = f.o[0] + i * f.dx; lat = f.o[1] + j * f.dy;
            PJ_COORD c = proj_trans(pj, PJ_FWD, proj_coord(lat, lon, 0, 0)); E = c.xy.y; N = c.xy.x;
        } else {
            E = f.o[0] + i * f.dx; N = f.o[1] + j * f.dy;
            PJ_COORD c = proj_trans(pj, PJ_INV, proj_coord(N, E, 0, 0)); lat = c.lp.lam; lon = c.lp.phi;
        }
        const double z = f.o[2] + k * f.dz;
        double s[3] = {RX, RY, RZ}, e[3] = {E + offx, N + offy, z};
        const bool occ = T->Intersect3d(s, e, IntersectCallback);
        const float want = occ ? 0.f : (vfloat ? KR(lon, lat, z) : (float)KC(lon, lat, z));
        const float got = f.at(i, j, k);
        nvis += !occ;
        if (memcmp(&want, &got, 4) != 0) {
            if (mism < 10) printf("  mismatch at (%d,%d,%d): recomputed %.9g, file %.9g\n", i, j, k, want, got);
            ++mism;
        }
    }
    proj_destroy(pj); proj_context_destroy(ctx);
    printf("random nodes %zu (visible %zu): bit-exact mismatches %zu\n", n, nvis, mism);
    return mism == 0 ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: field_check stats PREFIX | csv PREFIX CSV [TOL] | random PREFIX INDEX [N] [SEED]\n");
        return 2;
    }
    const std::string mode = argv[1], prefix = argv[2];
    try {
        if (mode == "stats") return stats(prefix);
        if (mode == "csv" && argc >= 4) return csv(prefix, argv[3], argc > 4 ? atof(argv[4]) : 1.2e-10);
        if (mode == "random" && argc >= 4)
            return random_check(prefix, argv[3], argc > 4 ? (size_t)atoll(argv[4]) : 10000, argc > 5 ? strtoull(argv[5], 0, 10) : 1);
    } catch (const std::exception& e) {
        fprintf(stderr, "error: %s\n", e.what()); return 1;
    }
    fprintf(stderr, "bad arguments\n");
    return 2;
}
