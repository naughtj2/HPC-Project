CXX ?= g++
CXXFLAGS ?= -O3 -std=gnu++17 -pthread
LDFLAGS ?= -pthread


SRC := $(wildcard *.cc)
OBJ := $(SRC:.cc=.o)

all: pricer

pricer: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ) pricer
