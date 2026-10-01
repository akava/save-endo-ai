CXX?=clang++
CXXFLAGS=-O2 -std=c++17
all: build/dna build/build
build/dna: src/dna.cpp
	@mkdir -p build out
	$(CXX) $(CXXFLAGS) -o $@ $<
build/build: src/build.cpp
	@mkdir -p build out
	$(CXX) $(CXXFLAGS) -o $@ $< -lz
data/endo.dna: doc/endo.zip
	unzip -o -q doc/endo.zip -d data
