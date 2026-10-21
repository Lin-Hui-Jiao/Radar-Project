// site_check.cpp — checks for the vehicle-radar 512^3 fields, using the project's R-tree code (read only):
//
//   site_check align  IDX_A IDX_B OFFX OFFY E0 N0 NX NY
//       Top-surface heights (Getheight3d) of two indexes on a 1 m EPSG:2326 grid; where both have triangles they must
//       be identical if both indexes share the R-tree frame local = EPSG:2326 + (OFFX, OFFY).
//   site_check cover  IDX OFFX OFFY E0 N0 NX NY
//       Columns of the window without any triangle (no ground, no buildings) and their extent.
//   site_check heights IDX OFFX OFFY E0 N0 NX NY OUT.f32
//       Writes the top-surface map (float32 [N][E], -1e9 where no triangle), e.g. for the slice figures.
//   site_check section IDX OFFX OFFY ew|ns|h C A0 A1 B0 B1 OUT.txt
//       True mesh cross-section of a plane: one line "a0 b0 a1 b1" per triangle cut by the plane (EPSG:2326 metres
//       and absolute altitude). ew: plane N = C, a = E in [A0, A1], b = z in [B0, B1]; ns: plane E = C, a = N, b = z;
//       h: plane z = C, a = E in [A0, A1], b = N in [B0, B1]. Unlike the top surface it shows overhangs and recesses.
//   site_check site   IDX OFFX OFFY E N MOUNT [SEARCH_M]
//       Is (E, N) a place where a vehicle can stand? Rules (all must hold):
//         ground  = Getheight3d at (E, N) (highest surface there)
//         street  : ground <= min(top surface within 15 m) + 1.0 m        (not a roof, podium or flyover deck)
//         flat    : max - min of the top surface within 2 m <= 0.3 m      (vehicle footprint)
//         clear   : nothing between ground + 0.05 m and the antenna at ground + MOUNT
//         sky     : nothing above the antenna (vertical segment up to +1000 m), i.e. not under a bridge or a tree
//       With SEARCH_M > 0 the nodes within SEARCH_M metres (1 m steps) that pass are listed, nearest first.
// Exit code 0 = pass (align: identical; site: the given point passes), 1 = fail.
#include "hiradar/occlusion_utils.h"
#include <proj.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

bool IntersectCallback(ValueType id, const double a[3], const double b[3]);
double HeightCallback(ValueType id, const double x, const double y);

namespace {

std::unique_ptr<RTree3d> load(const char* f) {
    std::unique_ptr<RTree3d> T(new RTree3d());
    if (!T->Load(f)) { fprintf(stderr, "cannot load %s\n", f); exit(2); }
    return T;
}

// same axis convention as field_export (authority order: lat, lon in; northing in xy.x, easting in xy.y)
void to_lonlat(double E, double N, double& lon, double& lat) {
    PJ_CONTEXT* ctx = proj_context_create();
    PJ* pj = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:2326", NULL);
    PJ_COORD c = proj_trans(pj, PJ_INV, proj_coord(N, E, 0, 0));
    lat = c.lp.lam; lon = c.lp.phi;
    proj_destroy(pj); proj_context_destroy(ctx);
}

std::vector<double> heights(RTree3d& T, double offx, double offy, double E0, double N0, int nx, int ny) {
    std::vector<double> h((size_t)nx * ny);
    #pragma omp parallel for schedule(dynamic, 4)
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) h[(size_t)j * nx + i] = T.Getheight3d(E0 + i + offx, N0 + j + offy, HeightCallback);
    return h;
}

int align(int argc, char** argv) {
    if (argc != 10) return 2;
    const double offx = atof(argv[4]), offy = atof(argv[5]), E0 = atof(argv[6]), N0 = atof(argv[7]);
    const int nx = atoi(argv[8]), ny = atoi(argv[9]);
    std::vector<double> a, b;
    { auto T = load(argv[2]); a = heights(*T, offx, offy, E0, N0, nx, ny); }
    { auto T = load(argv[3]); b = heights(*T, offx, offy, E0, N0, nx, ny); }
    size_t both = 0, same = 0, d1 = 0;
    double dmax = 0;
    for (size_t t = 0; t < a.size(); ++t) {
        if (a[t] < MINH || b[t] < MINH) continue;
        ++both;
        const double d = std::fabs(a[t] - b[t]);
        same += d == 0; d1 += d > 1.0; dmax = std::max(dmax, d);
    }
    printf("align: offset (%.12f, %.12f), %zu nodes covered by both indexes: identical %zu, |diff| > 1 m %zu, max |diff| %.4f m\n",
           offx, offy, both, same, d1, dmax);
    return both > 0 && same == both ? 0 : 1;
}

int dump(int argc, char** argv) {   // top-surface map, float32 [N][E], -1e9 where no triangle
    if (argc != 10) return 2;
    const double offx = atof(argv[3]), offy = atof(argv[4]), E0 = atof(argv[5]), N0 = atof(argv[6]);
    const int nx = atoi(argv[7]), ny = atoi(argv[8]);
    auto T = load(argv[2]);
    const std::vector<double> h = heights(*T, offx, offy, E0, N0, nx, ny);
    std::vector<float> f(h.size());
    for (size_t t = 0; t < h.size(); ++t) f[t] = h[t] < MINH ? -1.0e9f : (float)h[t];
    FILE* o = fopen(argv[9], "wb");
    if (!o || fwrite(f.data(), 4, f.size(), o) != f.size()) { fprintf(stderr, "cannot write %s\n", argv[9]); return 1; }
    fclose(o);
    printf("heights: wrote %s (%d x %d, float32 [N][E], origin E %.0f N %.0f)\n", argv[9], nx, ny, E0, N0);
    return 0;
}

int section(int argc, char** argv) {
    if (argc != 12) return 2;
    const double off[3] = {atof(argv[3]), atof(argv[4]), 0.0};
    const std::string ax = argv[5];
    int p, u, v;   // plane axis and the two in-plane axes (0 = E, 1 = N, 2 = z)
    if (ax == "ew") { p = 1; u = 0; v = 2; }
    else if (ax == "ns") { p = 0; u = 1; v = 2; }
    else if (ax == "h") { p = 2; u = 0; v = 1; }
    else return 2;
    const double C = atof(argv[6]) + off[p];
    double mn[3], mx[3];
    mn[p] = C - 1e-6; mx[p] = C + 1e-6;
    mn[u] = atof(argv[7]) + off[u]; mx[u] = atof(argv[8]) + off[u];
    mn[v] = atof(argv[9]) + off[v]; mx[v] = atof(argv[10]) + off[v];
    auto T = load(argv[2]);
    FILE* o = fopen(argv[11], "w");
    if (!o) { fprintf(stderr, "cannot write %s\n", argv[11]); return 1; }
    size_t ntri = 0, nseg = 0;
    T->Search(mn, mx, [&](const ValueType& t) {
        const double w[3][3] = {{std::get<0>(t), std::get<1>(t), std::get<2>(t)},
                                {std::get<3>(t), std::get<4>(t), std::get<5>(t)},
                                {std::get<6>(t), std::get<7>(t), std::get<8>(t)}};
        ++ntri;
        double pa[3], pb[3];
        int n = 0;
        for (int a = 0; a < 3 && n < 2; ++a) {   // vertices on the plane and edges crossing it
            const int b = (a + 1) % 3;
            const double da = w[a][p] - C, db = w[b][p] - C;
            if (da == 0) { pa[n] = w[a][u]; pb[n] = w[a][v]; ++n; }
            if (n < 2 && ((da < 0 && db > 0) || (da > 0 && db < 0))) {
                const double s = da / (da - db);
                pa[n] = w[a][u] + s * (w[b][u] - w[a][u]); pb[n] = w[a][v] + s * (w[b][v] - w[a][v]); ++n;
            }
        }
        if (n == 2) {
            fprintf(o, "%.4f %.4f %.4f %.4f\n", pa[0] - off[u], pb[0] - off[v], pa[1] - off[u], pb[1] - off[v]);
            ++nseg;
        }
        return true;
    });
    fclose(o);
    printf("section %s %s: %zu triangles in the box, %zu segments -> %s\n", ax.c_str(), argv[6], ntri, nseg, argv[11]);
    return 0;
}

int cover(int argc, char** argv) {
    if (argc != 9) return 2;
    const double offx = atof(argv[3]), offy = atof(argv[4]), E0 = atof(argv[5]), N0 = atof(argv[6]);
    const int nx = atoi(argv[7]), ny = atoi(argv[8]);
    auto T = load(argv[2]);
    const std::vector<double> h = heights(*T, offx, offy, E0, N0, nx, ny);
    size_t miss = 0;
    int ci0 = nx, ci1 = -1, cj0 = ny, cj1 = -1;   // extent of the covered columns
    double hmin = 1e300, hmax = -1e300;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double v = h[(size_t)j * nx + i];
            if (v < MINH) { ++miss; continue; }
            ci0 = std::min(ci0, i); ci1 = std::max(ci1, i); cj0 = std::min(cj0, j); cj1 = std::max(cj1, j);
            hmin = std::min(hmin, v); hmax = std::max(hmax, v);
        }
    printf("cover: window E %.0f..%.0f N %.0f..%.0f: %zu of %zu columns without triangles (%.2f%%); covered columns span "
           "i %d..%d (E %.0f..%.0f), j %d..%d (N %.0f..%.0f); top surface %.2f..%.2f m\n",
           E0, E0 + nx - 1, N0, N0 + ny - 1, miss, (size_t)nx * ny, 100.0 * miss / ((double)nx * ny), ci0, ci1, E0 + ci0,
           E0 + ci1, cj0, cj1, N0 + cj0, N0 + cj1, hmin, hmax);
    return 0;
}

struct SiteResult { double ground = 0, street = 0, flat = 0; bool covered = false, ok_street = false, ok_flat = false,
                    ok_clear = false, ok_sky = false; bool pass() const { return covered && ok_street && ok_flat && ok_clear && ok_sky; } };

SiteResult check_site(RTree3d& T, double offx, double offy, double E, double N, double mount) {
    SiteResult r;
    const double x = E + offx, y = N + offy;
    r.ground = T.Getheight3d(x, y, HeightCallback);
    if (r.ground < MINH) return r;
    r.covered = true;
    double mn15 = r.ground, lo = r.ground, hi = r.ground;
    for (int dj = -30; dj <= 30; ++dj)          // 0.5 m sampling
        for (int di = -30; di <= 30; ++di) {
            const double d = 0.5 * std::sqrt((double)(di * di + dj * dj));
            if (d > 15.0) continue;
            const double v = T.Getheight3d(x + 0.5 * di, y + 0.5 * dj, HeightCallback);
            if (v < MINH) continue;
            mn15 = std::min(mn15, v);
            if (d <= 2.0) { lo = std::min(lo, v); hi = std::max(hi, v); }
        }
    r.street = mn15; r.flat = hi - lo;
    r.ok_street = r.ground <= mn15 + 1.0;
    r.ok_flat = r.flat <= 0.3;
    double g[3] = {x, y, r.ground + 0.05}, a[3] = {x, y, r.ground + mount}, s[3] = {x, y, r.ground + mount + 1000.0};
    r.ok_clear = !T.Intersect3d(g, a, IntersectCallback);
    r.ok_sky = !T.Intersect3d(a, s, IntersectCallback);
    return r;
}

int site(int argc, char** argv) {
    if (argc != 8 && argc != 9) return 2;
    const double offx = atof(argv[3]), offy = atof(argv[4]), E = atof(argv[5]), N = atof(argv[6]), mount = atof(argv[7]);
    const int search = argc == 9 ? atoi(argv[8]) : 0;
    auto T = load(argv[2]);
    auto show = [&](double e, double n, const SiteResult& r) {
        double lon, lat; to_lonlat(e, n, lon, lat);
        printf("site E %.3f N %.3f (lon %.8f lat %.8f): ground %.3f m, street level (min within 15 m) %.3f m, "
               "flatness within 2 m %.3f m, antenna %.3f m -> street %s, flat %s, clear %s, sky %s => %s\n",
               e, n, lon, lat, r.ground, r.street, r.flat, r.ground + mount, r.ok_street ? "ok" : "NO", r.ok_flat ? "ok" : "NO",
               r.ok_clear ? "ok" : "NO", r.ok_sky ? "ok" : "NO", r.pass() ? "PASS" : "FAIL");
    };
    const SiteResult r0 = check_site(*T, offx, offy, E, N, mount);
    if (!r0.covered) printf("site E %.3f N %.3f: no triangle (outside the index)\n", E, N);
    else show(E, N, r0);
    if (search > 0) {
        struct C { double d; int di, dj; SiteResult r; };
        std::vector<C> pass;
        for (int dj = -search; dj <= search; ++dj)
            for (int di = -search; di <= search; ++di) {
                const double d = std::sqrt((double)(di * di + dj * dj));
                if (d > search || (di == 0 && dj == 0)) continue;
                const SiteResult r = check_site(*T, offx, offy, E + di, N + dj, mount);
                if (r.pass()) pass.push_back({d, di, dj, r});
            }
        std::sort(pass.begin(), pass.end(), [](const C& a, const C& b) { return a.d < b.d; });
        printf("passing nodes within %d m: %zu (nearest first)\n", search, pass.size());
        for (size_t t = 0; t < pass.size() && t < 8; ++t) {
            printf("  d %.2f m (dE %+d, dN %+d): ", pass[t].d, pass[t].di, pass[t].dj);
            show(E + pass[t].di, N + pass[t].dj, pass[t].r);
        }
    }
    return r0.pass() ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string m = argc > 1 ? argv[1] : "";
    int rc = 2;
    if (m == "align") rc = align(argc, argv);
    else if (m == "cover") rc = cover(argc, argv);
    else if (m == "heights") rc = dump(argc, argv);
    else if (m == "site") rc = site(argc, argv);
    else if (m == "section") rc = section(argc, argv);
    if (rc == 2)
        fprintf(stderr, "usage: site_check align IDX_A IDX_B OFFX OFFY E0 N0 NX NY\n"
                        "       site_check cover IDX OFFX OFFY E0 N0 NX NY\n"
                        "       site_check heights IDX OFFX OFFY E0 N0 NX NY OUT.f32\n"
                        "       site_check section IDX OFFX OFFY ew|ns|h C A0 A1 B0 B1 OUT.txt\n"
                        "       site_check site  IDX OFFX OFFY E N MOUNT [SEARCH_M]\n");
    return rc;
}
