###################################################
#
# Makefile for DESPOT Fitting Programs
# Tobias Wood December 2012
#
###################################################

# Set up Paths
SRC_DIR = .
BUILD_DIR := Build
INSTALL_DIR := .
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/system/gcc/gcc-4.8.0/lib64
EIGEN  = /software/local/k1078535/src/eigen/

$(INSTALL_BIN)/:
	mkdir -p $(INSTALL_BIN)
$(INSTALL_INC)/:
	mkdir -p $(INSTALL_INC)
$(INSTALL_LIB)/:
	mkdir -p $(INSTALL_LIB)

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/g++
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -msse3 -mssse3 -msse4.1 -msse4.2 -pthread -Wfatal-errors $(DEBUG) -DHAVE_NRECON
LD_FLAGS = -std=c++11 -lstdc++ -O3 -pthread -L$(INSTALL_LIB)
DEBUG_FLAGS = -g
RELEASE_FLAGS = -O3
INCLUDE = -I$(SRC_DIR) -I$(INSTALL_INC) -I$(EIGEN)
LD_LIBS  = -lrecon -lNifti -lz

#
# Pattern to build a .c in SRC_DIR/subdir to BUILD_DIR/subdir
# Note - put $(BUILD_DIR) last otherwise it gets prepended to all files including space
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	test -d $(BUILD_DIR) || mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< -o $@

#
# Variables to specify source files
#
DESPOT1_FILES = despot1_main.cpp DESPOT.cpp
DESPOT1_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(DESPOT1_FILES))
DESPOT-HIFI_FILES = despot_hifi.cpp DESPOT.cpp
DESPOT-HIFI_DEPS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(DESPOT-HIFI_FILES))
DESPOT2_FILES = despot2_main.cpp DESPOT.cpp
DESPOT2_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(DESPOT2_FILES))
MCDESPOT_FILES = mcdespot_main.cpp DESPOT.cpp
MCDESPOT_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(MCDESPOT_FILES))
THRESHOLD_FILES = threshold_main.cpp
THRESHOLD_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(THRESHOLD_FILES))
PHASEMAP_FILES = phasemap_main.cpp
PHASEMAP_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(PHASEMAP_FILES))
AFI_FILES      = afi_main.cpp
AFI_DEPS       = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(AFI_FILES))
CROP_FILES     = crop_main.cpp
CROP_DEPS      = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(CROP_FILES))

#
# Targets
#
despot1 : $(DESPOT1_DEPS)
	$(CXX) $(DESPOT1_DEPS) -o $(BUILD_DIR)/despot1 $(LD_FLAGS) $(LD_LIBS)
despot-hifi : $(DESPOT-HIFI_DEPS)
	$(CXX) $(DESPOT-HIFI_DEPS) -o $(BUILD_DIR)/despot-hifi $(LD_FLAGS) $(LD_LIBS)
despot2 : $(DESPOT2_DEPS)
	$(CXX) $(DESPOT2_DEPS) -o $(BUILD_DIR)/despot2 $(LD_FLAGS) $(LD_LIBS)
mcdespot : $(MCDESPOT_DEPS)
	$(CXX) $(MCDESPOT_DEPS) -o $(BUILD_DIR)/mcdespot $(LD_FLAGS) $(LD_LIBS)
phasemap : $(PHASEMAP_DEPS)
	$(CXX) $(PHASEMAP_DEPS) -o $(BUILD_DIR)/phasemap $(LD_FLAGS) $(LD_LIBS)
afi : $(AFI_DEPS)
	$(CXX) $(AFI_DEPS) -o $(BUILD_DIR)/afi $(LD_FLAGS) $(LD_LIBS)

.PHONY : build release debug install clean

release : CXX_FLAGS+=$(RELEASE_FLAGS)
release : build

debug   : CXX_FLAGS+=$(DEBUG_FLAGS)
debug   : build

build : despot1 despot2 despot-hifi mcdespot phasemap afi

install : $(INSTALL_BIN)
	cp $(BUILD_DIR)/despot1 $(INSTALL_BIN)/
	cp $(BUILD_DIR)/despot2 $(INSTALL_BIN)/
	cp $(BUILD_DIR)/despot-hifi $(INSTALL_BIN)/
	cp $(BUILD_DIR)/mcdespot $(INSTALL_BIN)/
	cp $(BUILD_DIR)/phasemap $(INSTALL_BIN)/
	cp $(BUILD_DIR)/afi $(INSTALL_BIN)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
