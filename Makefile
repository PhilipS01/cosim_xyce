CXX := clang++
CXXFLAGS := -std=c++17 -Iinclude -Wall -Wextra

# Main partitioned solver (excludes the GenExt driver, which has its own main()).
SRC := $(filter-out src/GenExtDriver.cpp,$(wildcard src/*.cpp))
OUT := main

# Xyce General External Device driver (links libxyce).
XYCE_PREFIX := /usr/local/XyceNF_7.10
GENEXT_SRC  := src/GenExtDriver.cpp
GENEXT_OUT  := genext_driver

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

genext: $(GENEXT_SRC)
	$(CXX) -std=c++17 -I$(XYCE_PREFIX)/include $(GENEXT_SRC) \
		-L$(XYCE_PREFIX)/lib -lxyce -Wl,-rpath,$(XYCE_PREFIX)/lib \
		-framework Accelerate -o $(GENEXT_OUT)

run: all
	./$(OUT)

clean:
	rm -f $(OUT) $(GENEXT_OUT)
