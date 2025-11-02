echo clean up
rm dist/*
cd radar && make clean && cd -
make clean
echo sync files
cp radar_client.cpp radar/
cp -r vendor radar/
cp -r db radar/
cp include/hiradar/radar_pool.hpp radar/
echo distribute 8 radar
tar czf dist/8_radar.tar.gz -C radar .
echo distribute geosot3d
tar czf dist/geosot3d.tar.gz -C geosot3d .
echo distribute general radar
tar --exclude-vcs --exclude='radar' --exclude='geosot3d' --exclude='dist' -czf dist/general_radar.tar.gz .
