# Compiler settings
CXX = g++

PPSS_TOOLS_PATH = /home/p/Documents/PPSS_tools

CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -I./include -I$(PPSS_TOOLS_PATH) -I$(PPSS_TOOLS_PATH)/Tools -I$(PPSS_TOOLS_PATH)/Tools/src/include -I$(PPSS_TOOLS_PATH)/BetheBloch/kit-dedx-fitter/src -I$(PPSS_TOOLS_PATH)/BetheBloch/src/include

# ROOT settings
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS += $(ROOTCFLAGS)

# --- LINKING LIBRARIES ---
# Direct path to PPSS shared libraries
PPSS_LIBS = $(PPSS_TOOLS_PATH)/Tools/Event.so \
            $(PPSS_TOOLS_PATH)/Tools/CutsMap.noDict.so \
            $(PPSS_TOOLS_PATH)/Tools/EventMixer.noDict.so \
            $(PPSS_TOOLS_PATH)/Tools/EventXeLa.so \
            $(PPSS_TOOLS_PATH)/Tools/EventXeLaMag.so

# LDFLAGS contains tool libraries and Boost libraries
LDFLAGS = $(PPSS_LIBS) -lboost_context -lboost_coroutine
# ---------------------------------------

# Directories
SRC_DIR = src
CORE_DIR = src/core
DIR_40 = src/40
DIR_75 = src/75
INC_DIR = include
BUILD_DIR = build

# External source files to compile along
BETHE_BLOCH_SRC = $(PPSS_TOOLS_PATH)/BetheBloch/src/BetheBlochWrapper.cc

# 1. Common files
COMMON_SRCS = $(filter-out $(DIR_40)/main%.cpp $(DIR_40)/Analysis%.cpp $(SRC_DIR)/main%.cpp $(SRC_DIR)/Analysis%.cpp, \
              $(wildcard $(SRC_DIR)/*.cpp) \
              $(wildcard $(CORE_DIR)/*.cpp) \
              $(wildcard $(DIR_40)/*.cpp) \
              $(wildcard $(DIR_75)/*.cpp))

COMMON_OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(notdir $(COMMON_SRCS))) $(BUILD_DIR)/BetheBlochWrapper.o

# 2. Executable specific objects
DATA_OBJS = $(BUILD_DIR)/main.o $(BUILD_DIR)/Analysis.o
MIXED_OBJS = $(BUILD_DIR)/main_mixed.o $(BUILD_DIR)/Analysis_mixed.o
COMPARE_OBJS = $(BUILD_DIR)/main_compare.o
EVENTS_OBJS = $(BUILD_DIR)/main_events.o $(BUILD_DIR)/Analysis_Events.o

VPATH = $(SRC_DIR):$(CORE_DIR):$(DIR_40):$(DIR_75)

# 3. Targets
TARGET_DATA = analysis
TARGET_MIXED = run_mixed
TARGET_COMPARE = compare_f2
TARGET_EVENTS = run_events

all: $(BUILD_DIR) $(TARGET_DATA) $(TARGET_MIXED) $(TARGET_COMPARE) $(TARGET_EVENTS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Added $(LDFLAGS) at the end of linking commands
$(TARGET_DATA): $(COMMON_OBJS) $(DATA_OBJS)
	$(CXX) $^ $(ROOTLIBS) $(LDFLAGS) -o $@

$(TARGET_MIXED): $(COMMON_OBJS) $(MIXED_OBJS)
	$(CXX) $^ $(ROOTLIBS) $(LDFLAGS) -o $@

$(TARGET_COMPARE): $(COMMON_OBJS) $(COMPARE_OBJS)
	$(CXX) $^ $(ROOTLIBS) $(LDFLAGS) -o $@

$(TARGET_EVENTS): $(COMMON_OBJS) $(EVENTS_OBJS)
	$(CXX) $^ $(ROOTLIBS) $(LDFLAGS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/BetheBlochWrapper.o: $(BETHE_BLOCH_SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET_DATA) $(TARGET_MIXED) $(TARGET_COMPARE) $(TARGET_EVENTS)

.PHONY: all clean