CXX=g++
CXX_FLAGS=-O2 -g -DDEBUG -std=c++17 -I./include -I/usr/local/include -I/usr/include/gdal -L/usr/local/lib 
LIBS = -lgeosot3d -lclickhouse-cpp-lib -lpthread -llz4 -lzstd -lcityhash -lgdal
BUILD=build
PREFIX=/usr/local/

all:  $(BUILD)/libhiradar.a $(BUILD)/libhiradar.so

$(BUILD)/%.o: %.cpp
	mpic++ $(CXX_FLAGS) -fPIC -c -o $@ $^

$(BUILD)/libhiradar.a: $(BUILD)/caltools.o $(BUILD)/radar.o $(BUILD)/TerrainGrid.o
	$(AR) -cr $@ $^

$(BUILD)/libhiradar.so: $(BUILD)/caltools.o $(BUILD)/radar.o $(BUILD)/TerrainGrid.o
	mpic++ -shared -fPIC $(CXX_FLAGS) -o $@ $^

# $(BUILD)/radar_server: $(BUILD)/radar_server.o $(BUILD)/libhiradar.a
# 	mpic++ $(CXX_FLAGS) -o $@ $^ $(LIBS)

# $(BUILD)/radar_client: radar_client.cpp
# 	$(CXX) $(CXX_FLAGS) -o $@ $^ $(LIBS)

# $(BUILD)/main: $(BUILD)/main.o $(BUILD)/libhiradar.a 
# 	mpic++ $(CXX_FLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf build/*

.PHONY: install
install:
	cp $(BUILD)/libhiradar.* $(PREFIX)/lib/
	cp -r include/* $(PREFIX)/include/
