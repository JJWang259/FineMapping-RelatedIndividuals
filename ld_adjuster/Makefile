# Makefile for LD Adjuster

CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native -DNDEBUG -fopenmp
INCLUDES = -I/usr/include/eigen3
LDFLAGS = -fopenmp
LIBS = 

# Detect if OpenBLAS is available
ifneq ($(shell pkg-config --exists openblas && echo "yes"),)
    LIBS += $(shell pkg-config --libs openblas)
    CXXFLAGS += -DEIGEN_USE_BLAS
    $(info Using OpenBLAS)
else ifneq ($(shell pkg-config --exists blas && echo "yes"),)
    LIBS += $(shell pkg-config --libs blas)
    CXXFLAGS += -DEIGEN_USE_BLAS
    $(info Using system BLAS)
else
    $(info No optimized BLAS found, using Eigen's built-in routines)
endif

# Detect Eigen3 location
ifeq ($(shell test -d /usr/include/eigen3 && echo "yes"),yes)
    INCLUDES = -I/usr/include/eigen3
else ifeq ($(shell test -d /usr/local/include/eigen3 && echo "yes"),yes)
    INCLUDES = -I/usr/local/include/eigen3
else ifeq ($(shell pkg-config --exists eigen3 && echo "yes"),yes)
    INCLUDES = $(shell pkg-config --cflags eigen3)
else
    $(error Eigen3 not found. Please install eigen3-dev or libeigen3-dev)
endif

TARGET = ld_adjuster
SOURCE = ld_adjuster.cpp

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Debug version
debug: CXXFLAGS = -std=c++11 -g -O0 -DDEBUG
debug: $(TARGET)

# Test with small data
test: $(TARGET)
	@echo "Testing with example data (you need to provide test files)"
	./$(TARGET) --raw test.raw --grm test --h2 0.5 --out test_output

# Show build information
info:
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Includes: $(INCLUDES)"
	@echo "Libraries: $(LIBS)"