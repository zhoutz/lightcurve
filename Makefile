executables := sd1a sd1b
headers := $(wildcard src/*)

ifneq (command line,$(origin CXX))
  CXX := clang++
endif

CXXFLAGS := -std=c++20 -Isrc -O3 -march=native #-ffast-math

all: $(executables)

$(executables): % : src/%.cpp $(headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
