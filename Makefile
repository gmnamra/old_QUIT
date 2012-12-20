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
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/local/k1078535
EIGEN  = /software/local/k1078535/include/eigen3

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
LD_FLAGS = -std=c++11 -stdlib=libc++ -O3 -L$(LIBCPP)/lib
INCLUDE = -I$(SRC_DIR) -I$(EIGEN) -cxx-isystem$(LIBCPP)/include/c++/v1

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

install : all $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp $(SRC_DIR)/*.h $(INSTALL_INC)/
	cp $(BUILD_DIR)/libNiftiImage.a $(INSTALL_LIB)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
