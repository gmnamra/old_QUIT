###################################################
#
# Makefile for librecon
# Tobias Wood December 2012
#
###################################################

#
# Macros
#

# Set up Directories
SRC_DIR = .
BUILD_DIR := Build
INSTALL_DIR := $(BUILD_DIR)
LIB_DIR := $(INSTALL_DIR)/lib
INCLUDE_DIR := $(INSTALL_DIR)/include

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
$(LIB_DIR):
	mkdir -p $(LIB_DIR)
$(INCLUDE_DIR):
	mkdir -p $(INCLUDE_DIR)

# Set up all our compiler options
CXX = clang++
AR = ar rcs
CXX_FLAGS = -std=c++11 -stdlib=libc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2
LD_FLAGS = -std=c++11 -stdlib=libc++ -O3 -L$(BUILD_DIR) -L/home/k1078535/local/lib
INCLUDE = -I$(SRC_DIR) -I/home/k1078535/local/include/cxx

#
# Pattern to build a .c in SRC_DIR/subdir to BUILD_DIR/subdir
# Note - put $(BUILD_DIR) last otherwise it gets prepended to all files including space
$(BUILD)/%.o : $(SRC_DIR)/%.cpp $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< -o $@

#
# Variables to actually specify source files
#

RECON_FILES = recon_util.cpp procpar.cpp FIDFile.cpp FID.cpp
RECON_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(RECON_FILES))
PROCPARSE_FILES = main_procparse.cpp
PROCPARSE_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(PROCPARSE_FILES))

all     : librecon procparse

install : all
	cp $(SRC_DIR)/*.h $(INCLUDE_DIR)
	cp $(BUILD_DIR)/librecon.a $(LIB_DIR)
	cp $(BUILD_DIR)/procparse $(INSTALL_DIR)

librecon : $(RECON_DEPS)
	$(AR) $(LIB_OPTIONS) $(BUILD_DIR)/librecon.a $(RECON_DEPS)

procparse : $(PROCPARSE_DEPS) librecon
	$(CXX) $(PROCPARSE_DEPS) -o $(BUILD_DIR)/procparse $(LD_FLAGS) -lrecon

clean : 
	rm -rf $(BUILD_DIR)/*
# For debugging variables
print-%:
	@echo $* = $($*)
