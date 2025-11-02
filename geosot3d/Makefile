PLATFORM=$(shell uname)
BUILD=./build
CXX_FLAGS=-O2 -g -DDEBUG -std=c++17 -fPIC
PREFIX=/usr/local
SHARED.Darwin=$(BUILD)/libgeosot3d.dylib
SHARED.Linux=$(BUILD)/libgeosot3d.so
SHARED=$(SHARED.$(PLATFORM))
LIBS=


all: $(BUILD)/geosot3d

$(BUILD)/%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $^

$(BUILD)/geosot3d: $(BUILD)/demo.o $(BUILD)/geosot3d.o
	$(CXX) -o $@ $^ $(LIBS)

$(SHARED): $(BUILD)/geosot3d.o
	$(CXX) -shared -fPIC -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf build/*

.PHONY: install
install: $(SHARED)
	cp $(SHARED) $(PREFIX)/lib/
	cp geosot3d.hpp $(PREFIX)/include/

