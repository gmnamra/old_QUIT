###################################################
#
# Makefile for librecon
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

$(INSTALL_BIN)/:
	mkdir -p $(INSTALL_BIN)
$(INSTALL_INC)/:
	mkdir -p $(INSTALL_INC)
$(INSTALL_LIB)/:
	mkdir -p $(INSTALL_LIB)

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/g++
AR = LD_RUN_PATH=$(LIBCPP) ar rcs
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
LD_FLAGS = -std=c++11 -lstdc++ -O3 -L$(BUILD_DIR) -L${LIBCPP}
INCLUDE = -I$(SRC_DIR)

#
# Pattern to build a .c in SRC_DIR/subdir to BUILD_DIR/subdir
# Note - put $(BUILD_DIR) last otherwise it gets prepended to all files including space
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	test -d $(BUILD_DIR) || mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< -o $@

#
# Variables to specify source files
#

RECON_FILES = recon_util.cpp procpar.cpp FIDFile.cpp FID.cpp
RECON_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(RECON_FILES))
PROCPARSE_FILES = main_procparse.cpp
PROCPARSE_DEPS  = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(PROCPARSE_FILES))

#
# Targets
#

librecon : $(RECON_DEPS)
	$(AR) $(BUILD_DIR)/librecon.a $(RECON_DEPS)

procparse : $(PROCPARSE_DEPS) librecon
	$(CXX) $(PROCPARSE_DEPS) -o $(BUILD_DIR)/procparse $(LD_FLAGS) -lrecon

all     : librecon procparse

install : all $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp $(SRC_DIR)/*.h $(INSTALL_INC)/
	cp $(BUILD_DIR)/librecon.a $(INSTALL_LIB)/
	cp $(BUILD_DIR)/procparse $(INSTALL_BIN)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
