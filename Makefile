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
INC_DIR = /software/local/k1078535/include 
LIBCPP  = /software/local/k1078535
EIGEN   = /software/local/k1078535/include/eigen3

$(INSTALL_BIN)/:
	mkdir -p $(INSTALL_BIN)
$(INSTALL_INC)/:
	mkdir -p $(INCLUDE_DIR)
$(INSTALL_LIB)/:
	mkdir -p $(INSTALL_LIB)

# Set up all our compiler options
CXX = clang++
AR = ar rcs
CXX_FLAGS = -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -std=c++11 -stdlib=libc++ $(DEBUG)
LD_FLAGS = -std=c++11 -stdlib=libc++ -O3 -L$(INSTALL_LIB) -L$(LIBCPP)/lib
INCLUDE = -I$(SRC_DIR) -I$(INC_DIR) -I$(EIGEN) -cxx-isystem$(LIBCPP)/include/c++/v1
LD_LIBS  = -lrecon -lNiftiImage -lz

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

all     : despot1 despot2 despot-hifi mcdespot phasemap afi

install : all $(INSTALL_BIN)
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
