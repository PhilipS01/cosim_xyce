CXX := clang++
# -O2: the driver resamples every waveform once per WR iteration (windows x WR steps),
# which is a real share of the runtime; the build was unoptimized until now.
CXXFLAGS := -std=c++17 -O2 -Iinclude -Wall -Wextra

SRC := $(wildcard src/*.cpp)
OUT := main

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

run: all
	./$(OUT)

clean:
	rm -f $(OUT)
