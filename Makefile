v###################################################
#
# Makefile for NiftiImage
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
LIBCPP = /software/system/gcc/gcc-4.8.0/lib
EIGEN  = /software/local/k1078535/include/eigen3

$(INSTALL_BIN)/:
	mkdir -p $(INSTALL_BIN)
$(INSTALL_INC)/:
	mkdir -p $(INSTALL_INC)
$(INSTALL_LIB)/:
	mkdir -p $(INSTALL_LIB)

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/gcc
AR = LD_RUN_PATH=$(LIBCPP) ar rcs
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
LD_FLAGS = -O3 -L$(BUILD_DIR)
INCLUDE = -I$(SRC_DIR) -I$(EIGEN)

#
# Pattern to build a .c in SRC_DIR/subdir to BUILD_DIR/subdir
# Note - put $(BUILD_DIR) last otherwise it gets prepended to all files including space
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	test -d $(BUILD_DIR) || mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< -o $@

#
# Variables to specify source files
#

NIFTI_FILES = Nifti.cpp
NIFTI_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(NIFTI_FILES))
HDR_FILES   = nifti_hdr.cpp
HDR_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(HDR_FILES))

#
# Targets
#

Nifti : $(NIFTI_DEPS)
	$(AR) $(LIB_OPTIONS) $(BUILD_DIR)/libNifti.a $(NIFTI_DEPS)

nifti_hdr : $(HDR_DEPS) Nifti
	$(CXX) $(CXX_FLAGS) $(HDR_DEPS) -o $(BUILD_DIR)/nifti_hdr $(LD_FLAGS) -lNifti -lz

all     : Nifti nifti_hdr

install : all $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp $(SRC_DIR)/*.h $(INSTALL_INC)/
	cp $(BUILD_DIR)/libNifti.a $(INSTALL_LIB)/
	cp $(BUILD_DIR)/nifti_hdr $(INSTALL_BIN)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
