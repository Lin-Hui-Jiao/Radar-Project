PLATFORM=$(shell uname)
INCLUDE.Darwin=-I${JAVA_HOME}/include -I${JAVA_HOME}/include/darwin
INCLUDE.LInux=-I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux
INCLUDE=$(INCLUDE.$(PLATFORM))
BUILD=build
CXX=g++
CXX_FLAGS=-O2 -g -DDEBUG -std=c++17 -I./include -fPIC
PACKAGE_DIR=cn/org/higis/radar/
LIBS = -lgeosot3dr
SRC=$(addprefix $(PACKAGE_DIR),GeoSOT3D.java)

SHARED.Darwin=$(BUILD)/libgeosot3djni.dylib
SHARED.Linux=$(BUILD)/libgeosot3djni.so
SHARED=$(SHARED.$(PLATFORM))

$(shell mkdir -p $(BUILD) $(PACKAGE_DIR))

$(BUILD)/%.o: %.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $^

$(SHARED): $(BUILD)/cn_org_higis_radar_GeoSOT3D.o
	$(CXX) -o $@ -fPIC -shared -o $@ $^ -lc $(LIBS)

.PHONY: demo
demo: $(PACKAGE_DIR)/Demo.java $(SRC) $(SHARED)
	javac -d $(BUILD) $(PACKAGE_DIR)/Demo.java $(SRC)
	cd $(BUILD) && java -Djava.library.path=./ cn.org.higis.radar.Demo && cd -

.PHONY: gen_header
gen_header:
	javac -d $(BUILD) -h . $(SRC)

.PHONY: clean
clean:
	rm -rf build/*

