executables := sd1a sd1b sd1c sd1d sd1e sd1f sd1bdef
headers := $(wildcard src/*)

ifneq (command line,$(origin CXX))
  CXX := clang++
endif

CXXFLAGS := -std=c++20 -Isrc -O2 -march=native -ffast-math

all: $(executables)

$(executables): % : src/%.cpp $(headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
