# macOS ships clang++, MSYS2/MinGW ships g++ -- pick whichever exists so that a bare
# `make` works on both. Override explicitly with `make CXX=...` (setup.sh/setup.ps1 do).
CXX := $(shell if command -v clang++ >/dev/null 2>&1; then echo clang++; else echo g++; fi)
# -O2: the driver resamples every waveform once per WR iteration (windows x WR steps),
# which is a real share of the runtime; the build was unoptimized until now.
CXXFLAGS := -std=c++17 -O2 -Iinclude -Wall -Wextra

SRC := $(wildcard src/*.cpp)
OUT := main

# Windows: MSYS2's g++ links against libgcc_s_seh-1.dll / libstdc++-6.dll /
# libwinpthread-1.dll, which live in ucrt64\bin and are only on PATH inside an MSYS2
# shell (setup.ps1 adds them for its own session only). A dynamically linked main.exe
# then dies with "libgcc_s_seh-1.dll was not found" as soon as anything else launches
# it -- sim_ui.py included. Link the runtime in so the binary stands alone.
ifeq ($(OS),Windows_NT)
LDFLAGS := -static
endif

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT) $(LDFLAGS)

run: all
	./$(OUT)

clean:
	rm -f $(OUT)
