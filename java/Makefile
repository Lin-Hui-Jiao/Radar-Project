PLATFORM=$(shell uname)
INCLUDE.Darwin=-I${JAVA_HOME}/include -I${JAVA_HOME}/include/darwin
INCLUDE.LInux=-I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux
INCLUDE=$(INCLUDE.$(PLATFORM))
BUILD=build
CXX=mpic++
CXX_FLAGS=-O2 -g -DDEBUG -std=c++17 -I./include -I/usr/local/include -fPIC
PACKAGE_DIR=cn/org/higis/radar/
LIBS =-lgeosot3d -lclickhouse-cpp-lib -lpthread -lhiradar
SRC=$(addprefix $(PACKAGE_DIR),Radar.java Position.java)

SHARED.Darwin=$(BUILD)/libradarjni.dylib
SHARED.Linux=$(BUILD)/libradarjni.so
SHARED=$(SHARED.$(PLATFORM))

$(shell mkdir -p $(BUILD) $(PACKAGE_DIR))

$(BUILD)/%.o: %.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $^

$(SHARED): $(BUILD)/cn_org_higis_radar_Radar.o
	$(CXX) -o $@ -fPIC -shared -o $@ $^ -L/usr/local/lib -lc $(LIBS)

.PHONY: demo
demo: $(PACKAGE_DIR)/Demo.java $(SRC) $(SHARED)
	javac -d $(BUILD) $(PACKAGE_DIR)/Demo.java $(SRC)
	cd $(BUILD) && mpiexec -n 8 java -Djava.library.path=./ cn.org.higis.radar.Demo && cd -

.PHONY: gen_header
gen_header:
	javac -d $(BUILD) -h . $(SRC)

.PHONY: clean
clean:
	rm -rf build/*

