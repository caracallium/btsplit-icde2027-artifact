CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
LDFLAGS ?=

ifeq ($(OS),Windows_NT)
EXE = .exe
STACK_FLAGS = -Wl,--stack,67108864
DYNAMIC_LIBS = -lpsapi
endif

.PHONY: all
all: btsplit$(EXE) btsplit_dynamic$(EXE)

btsplit$(EXE): src/btsplit.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(STACK_FLAGS)

btsplit_dynamic$(EXE): src/btsplit_dynamic.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(STACK_FLAGS) $(DYNAMIC_LIBS)
