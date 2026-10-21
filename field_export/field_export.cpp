// field_export.cpp — exports one radar's direct-path power-density field on a regular 3D grid in the PSHC hand-off
// format: <out>.f32 (little-endian float32, [z][y][x], x fastest, W/m^2, occluded = 0) and <out>.json.
//
// Occlusion: the project's R-tree + segment/triangle test (RTree3d::Intersect3d with IntersectCallback,
// src/occlusion_utils.cpp), radar absolute altitude = Getheight3d at the radar + mounting height (as Radar::Move).
// Amplitude: ReplayKernel (V-float, bit-identical to the float32 chain of Radar::CalculateSinglePointPowerDensity)
// or CleanKernel (V-double), see radar_kernels.hpp.
// Derived from pshc/prototypes/rise_review/gen_mesh3d.cpp (lon/lat grid) and gen_sample1024c.cpp / gen_sample1024d.cpp
// (EPSG:2326 metric grid). Existing files of this repository are not modified.
#include "hiradar/occlusion_utils.h"
#include "radar_kernels.hpp"
#include "sha256.hpp"
#include <omp.h>
#include <chrono>
#include <ctime>
#include <memory>

bool IntersectCallback(ValueType id, const double a[3], const double b[3]);
double HeightCallback(ValueType id, const double x, const double y);

namespace {

struct Args {
    std::string index, grid = "metric", version = "v-float", out;
    double ox = 0, oy = 0;             // grid origin: lon/lat (lonlat) or EPSG:2326 easting/northing (metric)
    int nx = 0, ny = 0, nz = 0;
    double sx = 1, sy = 1, sz = 1;     // spacing in metres (lonlat: converted to degrees with lat_ref)
    double z0 = 0;                     // altitude of layer k = 0 (m, absolute)
    bool has_lat_ref = false; double lat_ref = 0;
    double rx = 0, ry = 0, rh = 80;    // radar position (lon/lat if |rx| <= 360, else EPSG:2326) and height above ground
    double elev_min = -90, elev_max = 90;
    double off_x = -836974.6354, off_y = -815456.5722;   // R-tree local = EPSG:2326 + offset
    int threads = 0;
    bool ok = false;
};

void usage() {
    fprintf(stderr,
        "usage: field_export --index FILE.3idx --grid lonlat|metric --origin X Y --n NX NY NZ [--spacing SX SY SZ]\n"
        "                    --z0 Z --radar X Y [--radar-height H] [--elev-min DEG] [--elev-max DEG]\n"
        "                    [--version v-float|v-double] [--rtree-offset OX OY] [--lat-ref LAT|LAT0,LAT1]\n"
        "                    --out PREFIX [--threads N]\n"
        "  --grid lonlat : origin is lon/lat (deg); spacing in metres is converted with 111320 m/deg and\n"
        "                  cos(lat_ref) like FixedHeightComparator::RunExperiment (lat_ref default: grid centre;\n"
        "                  'A,B' means (A+B)/2)\n"
        "  --grid metric : origin is EPSG:2326 easting/northing (m)\n"
        "  --radar       : lon/lat if |X| <= 360, else EPSG:2326; height is above the R-tree ground\n"
        "  defaults      : spacing 1 1 1, radar-height 80, elevation -90..90, v-float,\n"
        "                  rtree-offset -836974.6354 -815456.5722 (master Radar::Move: -836975 -815457)\n");
}

Args parse(int argc, char** argv) {
    Args a;
    auto need = [&](int i, int n) { if (i + n >= argc) { fprintf(stderr, "missing value for %s\n", argv[i]); exit(2); } };
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        if (k == "--index") { need(i, 1); a.index = argv[++i]; }
        else if (k == "--grid") { need(i, 1); a.grid = argv[++i]; }
        else if (k == "--origin") { need(i, 2); a.ox = atof(argv[++i]); a.oy = atof(argv[++i]); }
        else if (k == "--n") { need(i, 3); a.nx = atoi(argv[++i]); a.ny = atoi(argv[++i]); a.nz = atoi(argv[++i]); }
        else if (k == "--spacing") { need(i, 3); a.sx = atof(argv[++i]); a.sy = atof(argv[++i]); a.sz = atof(argv[++i]); }
        else if (k == "--z0") { need(i, 1); a.z0 = atof(argv[++i]); }
        else if (k == "--lat-ref") {
            need(i, 1); std::string v = argv[++i]; size_t c = v.find(',');
            a.lat_ref = c == std::string::npos ? atof(v.c_str()) : (atof(v.substr(0, c).c_str()) + atof(v.substr(c + 1).c_str())) / 2.0;
            a.has_lat_ref = true;
        }
        else if (k == "--radar") { need(i, 2); a.rx = atof(argv[++i]); a.ry = atof(argv[++i]); }
        else if (k == "--radar-height") { need(i, 1); a.rh = atof(argv[++i]); }
        else if (k == "--elev-min") { need(i, 1); a.elev_min = atof(argv[++i]); }
        else if (k == "--elev-max") { need(i, 1); a.elev_max = atof(argv[++i]); }
        else if (k == "--version") { need(i, 1); a.version = argv[++i]; }
        else if (k == "--rtree-offset") { need(i, 2); a.off_x = atof(argv[++i]); a.off_y = atof(argv[++i]); }
        else if (k == "--out") { need(i, 1); a.out = argv[++i]; }
        else if (k == "--threads") { need(i, 1); a.threads = atoi(argv[++i]); }
        else if (k == "-h" || k == "--help") { usage(); exit(0); }
        else { fprintf(stderr, "unknown option %s\n", argv[i]); usage(); exit(2); }
    }
    a.ok = !a.index.empty() && !a.out.empty() && a.nx > 0 && a.ny > 0 && a.nz > 0 &&
           (a.grid == "lonlat" || a.grid == "metric") && (a.version == "v-float" || a.version == "v-double");
    return a;
}

// EPSG:4326 <-> EPSG:2326 with the axis order used by Radar::Move / gen_sample1024c (authority order: lat, lon in;
// northing in xy.x, easting in xy.y out). One instance per thread: PROJ objects are not thread-safe.
struct Proj {
    PJ_CONTEXT* ctx = nullptr; PJ* pj = nullptr;
    Proj() {
        ctx = proj_context_create(); pj = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:2326", NULL);
        if (!pj) { fprintf(stderr, "PROJ: cannot create EPSG:4326 -> EPSG:2326\n"); exit(1); }
    }
    ~Proj() { if (pj) proj_destroy(pj); if (ctx) proj_context_destroy(ctx); }
    void fwd(double lon, double lat, double& E, double& N) const {
        PJ_COORD c = proj_trans(pj, PJ_FWD, proj_coord(lat, lon, 0, 0)); E = c.xy.y; N = c.xy.x;
    }
    void inv(double E, double N, double& lon, double& lat) const {
        PJ_COORD c = proj_trans(pj, PJ_INV, proj_coord(N, E, 0, 0)); lat = c.lp.lam; lon = c.lp.phi;
    }
};

std::string shell_line(const std::string& cmd) {
    std::string s; FILE* p = popen(cmd.c_str(), "r"); if (!p) return s;
    char buf[256]; while (fgets(buf, sizeof buf, p)) s += buf;
    pclose(p);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

std::string json_str(const std::string& s) {
    std::string o = "\"";
    for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
    return o + "\"";
}

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

}  // namespace

int main(int argc, char** argv) {
    Args A = parse(argc, argv);
    if (!A.ok) { usage(); return 2; }
    if (A.threads > 0) omp_set_num_threads(A.threads);
    const int nthreads = omp_get_max_threads();
    const bool lonlat = A.grid == "lonlat", vfloat = A.version == "v-float";
    const bool radar_lonlat = std::fabs(A.rx) <= 360.0;

    auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<RTree3d> T(new RTree3d());
    if (!T->Load(A.index.c_str())) { fprintf(stderr, "cannot load R-tree index %s\n", A.index.c_str()); return 1; }
    const double t_load = seconds_since(t0);

    Proj P0;
    // ---- radar position: occlusion frame (R-tree local) and amplitude frame (lon/lat)
    double r_lon, r_lat, r_E, r_N;
    if (radar_lonlat) {
        r_lon = A.rx; r_lat = A.ry;
        double olon = r_lon, olat = r_lat;
        if (vfloat) {   // Radar::Move: lon/lat held as float degrees-minutes-seconds, DMSToDecimal in float
            float rl[3], rb[3];
            ReplayKernel::DecimalToDMS(rl, (float)r_lon); ReplayKernel::DecimalToDMS(rb, (float)r_lat);
            olon = rl[0] + rl[1] / 60.0f + rl[2] / 3600.0f; olat = rb[0] + rb[1] / 60.0f + rb[2] / 3600.0f;
        }
        P0.fwd(olon, olat, r_E, r_N);
    } else {
        r_E = A.rx; r_N = A.ry; P0.inv(r_E, r_N, r_lon, r_lat);
    }
    const double RX = r_E + A.off_x, RY = r_N + A.off_y;
    const double ground = T->Getheight3d(RX, RY, HeightCallback);
    if (ground < MINH) { fprintf(stderr, "radar (%.3f, %.3f) is outside the index (no ground triangle)\n", RX, RY); return 1; }
    const double RZ = ground + A.rh;

    Meta km; km.rlon = r_lon; km.rlat = r_lat; km.ralt = vfloat ? (double)(float)RZ : RZ;
    km.minphi = A.elev_min; km.maxphi = A.elev_max;
    ReplayKernel KR; CleanKernel KC;
    if (vfloat) KR.init(km); else KC.init(km);

    // ---- grid
    double dlon = 0, dlat = 0, lat_ref = 0;
    if (lonlat) {
        dlat = A.sy / 111320.0;
        lat_ref = A.has_lat_ref ? A.lat_ref : A.oy + (A.ny - 1) * dlat / 2.0;
        dlon = A.sx / (111320.0 * std::cos(lat_ref * M_PI / 180.0));   // FixedHeightComparator::MetersToDegrees
    }
    const size_t ncol = (size_t)A.nx * A.ny;
    std::vector<double> LX(ncol), LY(ncol), LON(ncol), LAT(ncol);
    auto t1 = std::chrono::steady_clock::now();
    #pragma omp parallel
    {
        Proj P;
        #pragma omp for schedule(static)
        for (long c = 0; c < (long)ncol; ++c) {
            int i = (int)(c % A.nx), j = (int)(c / A.nx);
            double E, N, lon, lat;
            if (lonlat) { lon = A.ox + i * dlon; lat = A.oy + j * dlat; P.fwd(lon, lat, E, N); }
            else { E = A.ox + i * A.sx; N = A.oy + j * A.sy; P.inv(E, N, lon, lat); }
            LX[c] = E + A.off_x; LY[c] = N + A.off_y; LON[c] = lon; LAT[c] = lat;
        }
    }
    const double t_proj = seconds_since(t1);
    printf("index %s loaded in %.1f s; radar given as %s: lon/lat (%.7f, %.7f), EPSG:2326 (%.3f, %.3f), R-tree local (%.4f, %.4f)\n",
           A.index.c_str(), t_load, radar_lonlat ? "lon/lat" : "EPSG:2326", r_lon, r_lat, r_E, r_N, RX, RY);
    printf("ground %.6f m + %.3f m = radar absolute altitude %.6f m (kernel %.9g); %s, elevation %.1f..%.1f deg\n",
           ground, A.rh, RZ, km.ralt, vfloat ? "V-float" : "V-double", A.elev_min, A.elev_max);

    // ---- field, computed in z-slabs (keeps memory bounded for 1024^3)
    const std::string f32 = A.out + ".f32";
    FILE* fo = fopen(f32.c_str(), "wb");
    if (!fo) { fprintf(stderr, "cannot write %s\n", f32.c_str()); return 1; }
    const int slab = std::max(1, std::min(A.nz, (int)((size_t)256 << 20) / (int)std::max<size_t>(1, ncol * 4)));
    std::vector<float> buf(ncol * slab);
    pshc::Sha256 sha;
    size_t nvisible = 0, nvis_zero = 0;
    auto t2 = std::chrono::steady_clock::now();
    for (int k0 = 0; k0 < A.nz; k0 += slab) {
        const int kn = std::min(slab, A.nz - k0);
        #pragma omp parallel for schedule(dynamic, 16) reduction(+ : nvisible, nvis_zero)
        for (long c = 0; c < (long)ncol; ++c) {
            double s[3] = {RX, RY, RZ};
            for (int kk = 0; kk < kn; ++kk) {
                const double z = A.z0 + (k0 + kk) * A.sz;
                double e[3] = {LX[c], LY[c], z};
                const bool occ = T->Intersect3d(s, e, IntersectCallback);
                float v = 0.f;
                if (!occ) {
                    v = vfloat ? KR(LON[c], LAT[c], z) : (float)KC(LON[c], LAT[c], z);
                    ++nvisible; nvis_zero += v == 0.f;
                }
                buf[(size_t)kk * ncol + c] = v;
            }
        }
        const size_t bytes = (size_t)kn * ncol * 4;
        if (fwrite(buf.data(), 1, bytes, fo) != bytes) { fprintf(stderr, "write failed\n"); return 1; }
        sha.update(buf.data(), bytes);
    }
    fclose(fo);
    const double t_field = seconds_since(t2);
    const std::string digest = sha.hex();
    const size_t total = ncol * A.nz;
    printf("field %d x %d x %d computed in %.2f s (%d threads, projection %.2f s): visible %.4f, visible but zero amplitude %zu\n",
           A.nx, A.ny, A.nz, t_field, nthreads, t_proj, (double)nvisible / total, nvis_zero);

    // ---- provenance
#ifdef FIELD_EXPORT_SRC_DIR
    const std::string repo = std::string(FIELD_EXPORT_SRC_DIR) + "/..";
#else
    const std::string repo = ".";
#endif
    const std::string commit = shell_line("git -C '" + repo + "' rev-parse HEAD 2>/dev/null");
    const std::string branch = shell_line("git -C '" + repo + "' rev-parse --abbrev-ref HEAD 2>/dev/null");
    const bool dirty = !shell_line("git -C '" + repo + "' status --porcelain --untracked-files=no 2>/dev/null").empty();
    std::string cmd; for (int i = 0; i < argc; ++i) { if (i) cmd += ' '; cmd += argv[i]; }
    char when[64]; time_t now = time(nullptr); strftime(when, sizeof when, "%Y-%m-%dT%H:%M:%S%z", localtime(&now));

    const std::string js = A.out + ".json";
    FILE* fj = fopen(js.c_str(), "w");
    if (!fj) { fprintf(stderr, "cannot write %s\n", js.c_str()); return 1; }
    fprintf(fj, "{\n  \"grid\": {\n");
    fprintf(fj, "    \"nx\": %d, \"ny\": %d, \"nz\": %d,\n", A.nx, A.ny, A.nz);
    if (lonlat) {
        fprintf(fj, "    \"dx\": %.17g, \"dy\": %.17g, \"dz\": %.17g,\n", dlon, dlat, A.sz);
        fprintf(fj, "    \"crs\": \"EPSG:4326\", \"units\": [\"deg\", \"deg\", \"m\"],\n");
        fprintf(fj, "    \"lat_ref\": %.17g, \"spacing_m\": [%.17g, %.17g, %.17g],\n", lat_ref, A.sx, A.sy, A.sz);
    } else {
        fprintf(fj, "    \"dx\": %.17g, \"dy\": %.17g, \"dz\": %.17g,\n", A.sx, A.sy, A.sz);
        fprintf(fj, "    \"crs\": \"EPSG:2326\", \"units\": [\"m\", \"m\", \"m\"],\n");
    }
    fprintf(fj, "    \"origin\": [%.17g, %.17g, %.17g],\n", A.ox, A.oy, A.z0);
    fprintf(fj, "    \"vertical\": \"absolute altitude of layer k = origin[2] + k * dz\"\n  },\n");
    fprintf(fj, "  \"value\": {\"quantity\": \"direct-path power density\", \"unit\": \"W/m^2\", \"dtype\": \"float32\", "
                "\"endian\": \"little\", \"order\": \"zyx\", \"occluded\": 0},\n");
    fprintf(fj, "  \"sha256\": \"%s\",\n", digest.c_str());
    fprintf(fj, "  \"provenance\": {\n");
    fprintf(fj, "    \"generator\": \"field_export\", \"generator_repo\": %s, \"generator_branch\": %s,\n",
            json_str(repo).c_str(), json_str(branch).c_str());
    fprintf(fj, "    \"generator_commit\": %s, \"generator_tracked_files_modified\": %s,\n",
            json_str(commit).c_str(), dirty ? "true" : "false");
    fprintf(fj, "    \"command\": %s,\n", json_str(cmd).c_str());
    fprintf(fj, "    \"created\": \"%s\", \"threads\": %d,\n", when, nthreads);
    fprintf(fj, "    \"index\": %s, \"rtree_offset\": [%.17g, %.17g],\n", json_str(A.index).c_str(), A.off_x, A.off_y);
    fprintf(fj, "    \"radar\": {\n");
    fprintf(fj, "      \"input\": [%.17g, %.17g], \"input_crs\": \"%s\",\n", A.rx, A.ry, radar_lonlat ? "EPSG:4326" : "EPSG:2326");
    fprintf(fj, "      \"lonlat\": [%.17g, %.17g], \"epsg2326\": [%.17g, %.17g], \"rtree_local\": [%.17g, %.17g],\n",
            r_lon, r_lat, r_E, r_N, RX, RY);
    fprintf(fj, "      \"ground_m\": %.17g, \"mount_height_m\": %.17g, \"abs_alt_m\": %.17g, \"abs_alt_kernel_m\": %.17g,\n",
            ground, A.rh, RZ, km.ralt);
    fprintf(fj, "      \"type\": \"CityGuardRadar\", \"N\": %d, \"d_m\": %.17g, \"f_MHz\": %.17g, \"G_dB\": %.17g, \"Pt_kW\": %.17g,\n",
            km.N, km.d, km.f, km.G, km.Pt);
    fprintf(fj, "      \"azimuth_range_deg\": [%.17g, %.17g], \"elevation_range_deg\": [%.17g, %.17g]\n    },\n",
            km.mintheta, km.maxtheta, km.minphi, km.maxphi);
    fprintf(fj, "    \"data_version\": \"%s\", \"amplitude_kernel\": \"%s\",\n", vfloat ? "V-float" : "V-double",
            vfloat ? "ReplayKernel (float32 chain of Radar::CalculateSinglePointPowerDensity)" : "CleanKernel (double precision)");
    fprintf(fj, "    \"occlusion\": \"RTree3d::Intersect3d + IntersectCallback (src/occlusion_utils.cpp), target position in double\",\n");
    fprintf(fj, "    \"grid_type\": \"%s\",\n", A.grid.c_str());
    fprintf(fj, "    \"stats\": {\"visible_fraction\": %.17g, \"visible_zero_amplitude\": %zu},\n", (double)nvisible / total, nvis_zero);
    fprintf(fj, "    \"seconds\": {\"index_load\": %.3f, \"projection\": %.3f, \"field\": %.3f}\n", t_load, t_proj, t_field);
    fprintf(fj, "  }\n}\n");
    fclose(fj);
    printf("wrote %s (%.2f MB, sha256 %s) and %s\n", f32.c_str(), total * 4.0 / 1e6, digest.c_str(), js.c_str());
    return 0;
}
