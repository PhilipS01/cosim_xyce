CXX := clang++
CXXFLAGS := -std=c++17 -Iinclude -Wall -Wextra

SRC := $(wildcard src/*.cpp)
OUT := main

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

run: all
	./$(OUT)

clean:
	rm -f $(OUT)
