# Compiler settings
CXX = g++

PPSS_TOOLS_PATH = /home/p/Documents/PPSS_tools

CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -I./include -I$(PPSS_TOOLS_PATH) -I$(PPSS_TOOLS_PATH)/Tools -I$(PPSS_TOOLS_PATH)/Tools/src/include -I$(PPSS_TOOLS_PATH)/BetheBloch/kit-dedx-fitter/src -I$(PPSS_TOOLS_PATH)/BetheBloch/src/include

# ROOT settings
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS += $(ROOTCFLAGS)

# Directories
SRC_DIR = src
CORE_DIR = src/core
DIR_40 = src/40
DIR_75 = src/75
INC_DIR = include
BUILD_DIR = build

# External source files to compile along
BETHE_BLOCH_SRC = $(PPSS_TOOLS_PATH)/BetheBloch/src/BetheBlochWrapper.cc

# Source files (find all .cpp files in specified directories)
SRCS = $(wildcard $(SRC_DIR)/*.cpp) \
       $(wildcard $(CORE_DIR)/*.cpp) \
       $(wildcard $(DIR_40)/*.cpp) \
       $(wildcard $(DIR_75)/*.cpp)

# Object files (placed in folder, including external wrapper)
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(notdir $(SRCS))) $(BUILD_DIR)/BetheBlochWrapper.o

# Tell Make where to look for source files
VPATH = $(SRC_DIR):$(CORE_DIR):$(DIR_40):$(DIR_75)

# Main executable name
TARGET = analysis

# Default rule
all: $(BUILD_DIR) $(TARGET)

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Link object files to create executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(ROOTLIBS) -o $(TARGET)

# Universal rule to compile any .cpp file found via VPATH
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Specific rule to compile external BetheBlochWrapper.cc
$(BUILD_DIR)/BetheBlochWrapper.o: $(BETHE_BLOCH_SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean