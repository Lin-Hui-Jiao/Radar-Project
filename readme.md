2km*2km
lon = 114.1838991865
lat = 22.2738043465






./fixed_height_experiment --index ../../test_area_2_2.3idx --min-x 836000 --min-y 813750 --index-range-x -974.6354 --index-range-y -1706.5722 --radar-lon 114.183899 --radar-lat 22.273804 --radar-alt 80 --min-lon 114.174274 --max-lon 114.193680 --min-lat 22.262560 --max-lat 22.280623 --fixed-height 30



它做的事情是：随机生成指定范围内的 10万 个经纬高点，雷达位置固定；每个点执行一次真实遮挡查询；callback 中统计实际进入精确三角形相交判断的候选面片数；一旦精确相交命中就 return false 提前结束，不再统计后续面片。没有进行雷达能量计算。

我没有修改 CMakeLists.txt 或其他源码。独立编译已通过：

rootjiao@localhost:~/radar-demo-0927-shiyan$ cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./random_visibility_candidate_experiment \
  --samples 1000000 \
  --index ../../test_area_2_2.3idx \
  --summary-only \
  --warmup-runs 1
========================================
  Random Visibility Candidate Experiment
========================================

[1/5] Loading RTree index...
  Index file: ../../test_area_2_2.3idx

[2/5] Preparing coordinate transform and radar point...
  Radar lon/lat: 114.1838991865, 22.2738043465
  Radar RTree xyz: 17.3648744574, -461.5721331790, 123.0998151692
  Radar relative altitude: 80.0000000000 m

[3/5] Generating and projecting random points...
  Samples: 1000000
  Lon range: [114.1742740000, 114.1936800000]
  Lat range: [22.2625600000, 22.2806230000]
  Alt range: [0.0000000000, 100.0000000000] m
  Seed: 20260517
  OpenMP threads available: 8
  Query timing mode: batch chrono only
  Result mode: summary-only maximum throughput
  Intersection kernel: fast inline Moller-Trumbore

  Random generation time: 0.0153495450 s
  Projection precompute time: 0.0954875200 s
  Valid prepared queries: 1000000
  Projection failures: 0

[warmup] Running 1 untimed core query pass(es)...
  warmup 1: 4.0127718270 s, 249204.3014435767 queries/s

[4/5] Running parallel visibility queries with early-stop exact intersection...
  Valid queries: 1000000
  Projection failures: 0

[5/5] Writing results and summary...
  CSV disabled (--no-csv)

========== Summary ==========
Threads used/requested: 8 / 0
Random generation time: 0.015350 s
Projection precompute time: 0.095488 s
Core visibility-query time: 0.518531 s
Core visibility-query throughput: 1928526.657550 queries/s
Average core query time: 0.518531 us
CSV write time: 0.000000 s
End-to-end measured time: 0.629368 s
End-to-end throughput: 1588896.483761 queries/s
Visible count: 164569
Occluded count: 835431
Total candidate triangles tested exactly: 2148659
Average candidates per query: 2.148659
Min candidates per query: 0
Max candidates per query: 54



rootjiao@localhost:~/radar-demo-0927-shiyan/build/bin$ ./random_visibility_candidate_experiment \  --samples 1000000 \  --index ../../test_area_1_1.3idx \
  --min-lon 114.157292795 \
  --max-lon 114.166996872 \
  --min-lat 22.273848477 \
  --max-lat 22.282879994 \
  --radar-lon 114.1612 \
  --radar-lat 22.2788 \
  --radar-alt 80 \
  --min-alt 0 \
  --max-alt 100 \
  --min-x 834250 \
  --min-y 815000 \
  --index-range-x -2724.635391 \
  --index-range-y -456.572217 \
  --summary-only \
  --warmup-runs 1
========================================
  Random Visibility Candidate Experiment
========================================

[1/5] Loading RTree index...
  Index file: ../../test_area_1_1.3idx

[2/5] Preparing coordinate transform and radar point...
  Radar lon/lat: 114.1612000000, 22.2788000000
  Radar RTree xyz: -2322.0200993772, 91.6832387759, 94.6938482878
  Radar relative altitude: 80.0000000000 m

[3/5] Generating and projecting random points...
  Samples: 1000000
  Lon range: [114.1572927950, 114.1669968720]
  Lat range: [22.2738484770, 22.2828799940]
  Alt range: [0.0000000000, 100.0000000000] m
  Seed: 20260517
  OpenMP threads available: 8
  Query timing mode: batch chrono only
  Result mode: summary-only maximum throughput
  Intersection kernel: fast inline Moller-Trumbore

  Random generation time: 0.0152298440 s
  Projection precompute time: 0.0828604140 s
  Valid prepared queries: 1000000
  Projection failures: 0

[warmup] Running 1 untimed core query pass(es)...
  warmup 1: 0.3035693310 s, 3294140.4084064080 queries/s

[4/5] Running parallel visibility queries with early-stop exact intersection...
  Valid queries: 1000000
  Projection failures: 0

[5/5] Writing results and summary...
  CSV disabled (--no-csv)

========== Summary ==========
Threads used/requested: 8 / 0
Random generation time: 0.015230 s
Projection precompute time: 0.082860 s
Core visibility-query time: 0.316222 s
Core visibility-query throughput: 3162331.481253 queries/s
Average core query time: 0.316222 us
CSV write time: 0.000000 s
End-to-end measured time: 0.414313 s
End-to-end throughput: 2413636.222088 queries/s
Visible count: 140319
Occluded count: 859681
Total candidate triangles tested exactly: 2031719
Average candidates per query: 2.031719
Min candidates per query: 0
Max candidates per query: 48






rootjiao@localhost:~/radar-demo-0927-shiyan/build/bin$ ./random_visibility_candidate_experiment \
  --samples 1000000 \
  --index ../../test_area.3idx \
  --min-lon 114.164571 \
  --max-lon 114.169422 \
  --min-lat 22.278365 \
  --max-lat 22.28288 \
  --radar-lon 114.1670 \
  --radar-lat 22.2806 \
  --radar-alt 80 \
  --min-alt 0 \
  --max-alt 100 \
  --min-x 835000 \
  --min-y 815500 \
  --index-range-x -1975 \
  --index-range-y 43 \
  --summary-only \
  --warmup-runs 1
========================================
  Random Visibility Candidate Experiment
========================================

[1/5] Loading RTree index...
  Index file: ../../test_area.3idx

[2/5] Preparing coordinate transform and radar point...
  Radar lon/lat: 114.1670000000, 22.2806000000
  Radar RTree xyz: -1724.6286079863, 290.5265729887, 164.4593203777
  Radar relative altitude: 80.0000000000 m

[3/5] Generating and projecting random points...
  Samples: 1000000
  Lon range: [114.1645710000, 114.1694220000]
  Lat range: [22.2783650000, 22.2828800000]
  Alt range: [0.0000000000, 100.0000000000] m
  Seed: 20260517
  OpenMP threads available: 8
  Query timing mode: batch chrono only
  Result mode: summary-only maximum throughput
  Intersection kernel: fast inline Moller-Trumbore

  Random generation time: 0.0150413220 s
  Projection precompute time: 0.0852377470 s
  Valid prepared queries: 1000000
  Projection failures: 0

[warmup] Running 1 untimed core query pass(es)...
  warmup 1: 0.1556616440 s, 6424190.1492444733 queries/s

[4/5] Running parallel visibility queries with early-stop exact intersection...
  Valid queries: 1000000
  Projection failures: 0

[5/5] Writing results and summary...
  CSV disabled (--no-csv)

========== Summary ==========
Threads used/requested: 8 / 0
Random generation time: 0.015041 s
Projection precompute time: 0.085238 s
Core visibility-query time: 0.152427 s
Core visibility-query throughput: 6560530.018397 queries/s
Average core query time: 0.152427 us
CSV write time: 0.000000 s
End-to-end measured time: 0.252706 s
End-to-end throughput: 3957171.094046 queries/s
Visible count: 672022
Occluded count: 327978
Total candidate triangles tested exactly: 827389
Average candidates per query: 0.827389
Min candidates per query: 0
Max candidates per query: 44