BUILD=./build
CXX_FLAGS=-O2 -g -DDEBUG -std=c++17 -I../include -I/usr/include/gdal
LIBS = -lgeosot3d -lclickhouse-cpp-lib -lpthread -lhiradar -llz4 -lzstd -lcityhash -lgdal
all: $(BUILD)/an_apg-66 $(BUILD)/an_aps-115 $(BUILD)/an_apy-1 $(BUILD)/an_apy-2 $(BUILD)/an_sps-49 $(BUILD)/an_sps-49 $(BUILD)/spy-1d_v $(BUILD)/an_aps-145 $(BUILD)/an_mpq-65 $(BUILD)/radar_client

$(BUILD)/an_apg-66: server.cpp
	mpic++ -DRADAR=ANAPG66 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_aps-145: server.cpp
	mpic++ -DRADAR=ANAPS145 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_aps-115: server.cpp
	mpic++ -DRADAR=ANAPS115 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_apy-1: server.cpp
	mpic++ -DRADAR=ANAPY12 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_apy-2: server.cpp
	mpic++ -DRADAR=ANAPY12 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_sps-49: server.cpp
	mpic++ -DRADAR=ANSPS49 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/spy-1d_v: server.cpp
	mpic++ -DRADAR=SPY1D $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/an_mpq-65: server.cpp
	mpic++ -DRADAR=ANMPQ65 $(CXX_FLAGS) -o $@ $^ $(LIBS)

$(BUILD)/radar_client: radar_client.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lpthread

.PHONY: clean
clean:
	rm -rf build/*

