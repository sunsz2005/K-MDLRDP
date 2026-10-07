CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Wpedantic
INCLUDES = -Iinclude -Iapps
COMMON = src/graph.cpp src/evaluator.cpp src/algorithms.cpp

.PHONY: all test clean

all: bin/kmdrdp_greedy.exe bin/kmdrdp_ga.exe bin/kmdrdp_vns.exe

bin/kmdrdp_greedy.exe: $(COMMON) apps/greedy_main.cpp
	@if not exist bin mkdir bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(COMMON) apps/greedy_main.cpp -o $@

bin/kmdrdp_ga.exe: $(COMMON) apps/ga_main.cpp
	@if not exist bin mkdir bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(COMMON) apps/ga_main.cpp -o $@

bin/kmdrdp_vns.exe: $(COMMON) apps/vns_main.cpp
	@if not exist bin mkdir bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(COMMON) apps/vns_main.cpp -o $@

bin/test_core.exe: $(COMMON) tests/test_core.cpp
	@if not exist bin mkdir bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(COMMON) tests/test_core.cpp -o $@

test: all bin/test_core.exe
	bin\test_core.exe
	python tests\test_exact.py
	python tests\test_hybrid.py

clean:
	@if exist bin rmdir /S /Q bin

