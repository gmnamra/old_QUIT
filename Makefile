###################################################
#
# Makefile for NiftiImage
# Tobias Wood December 2012
#
###################################################

# Set up Paths
SRC_DIR = .
BUILD_DIR := Build
INSTALL_DIR := .
BIN_DIR = $(INSTALL_DIR)/bin
INC_DIR = $(INSTALL_DIR)/include
LIB_DIR = $(INSTALL_DIR)/lib
LIBCPP = /software/local/k1078535
EIGEN  = /software/local/k1078535/include/eigen3

$(BIN_DIR)/:
	mkdir -p $(BIN_DIR)
$(INC_DIR)/:
	mkdir -p $(INCLUDE_DIR)
$(LIB_DIR)/:
	mkdir -p $(LIB_DIR)

# Set up all our compiler options
CXX = clang++
AR = ar rcs
CXX_FLAGS = -std=c++11 -stdlib=libc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2
LD_FLAGS = -std=c++11 -stdlib=libc++ -O3 -L$(BUILD_DIR) -L$(LIBCPP)/lib
INCLUDE = -I$(SRC_DIR) -I$(LIBCPP)/include/c++/v1 -I$(EIGEN)

#
# Pattern to build a .c in SRC_DIR/subdir to BUILD_DIR/subdir
# Note - put $(BUILD_DIR) last otherwise it gets prepended to all files including space
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	test -d $(BUILD_DIR) || mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< -o $@

#
# Variables to specify source files
#

NIFTI_FILES = NiftiImage.cpp
NIFTI_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(NIFTI_FILES))

#
# Targets
#

NiftiImage : $(NIFTI_DEPS)
	$(AR) $(LIB_OPTIONS) $(BUILD_DIR)/libNiftiImage.a $(NIFTI_DEPS)

all     : NiftiImage

install : all $(BIN_DIR) $(INC_DIR) $(LIB_DIR)
	cp $(SRC_DIR)/*.h $(INC_DIR)/
	cp $(BUILD_DIR)/libNiftiImage.a $(LIB_DIR)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
