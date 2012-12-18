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
BIN_DIR = $(INSTALL_DIR)/bin
INC_DIR = $(INSTALL_DIR)/include
LIB_DIR = $(INSTALL_DIR)/lib
LIBCPP = /software/local/k1078535

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
INCLUDE = -I$(SRC_DIR) -I$(LIBCPP)/include/c++/v1

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
	$(AR) $(LIB_OPTIONS) $(BUILD_DIR)/librecon.a $(RECON_DEPS)

procparse : $(PROCPARSE_DEPS) librecon
	$(CXX) $(PROCPARSE_DEPS) -o $(BUILD_DIR)/procparse $(LD_FLAGS) -lrecon

all     : librecon procparse

install : all $(BIN_DIR) $(INC_DIR) $(LIB_DIR)
	cp $(SRC_DIR)/*.h $(INC_DIR)/
	cp $(BUILD_DIR)/librecon.a $(LIB_DIR)/
	cp $(BUILD_DIR)/procparse $(BIN_DIR)/

clean : 
	rm -rf $(BUILD_DIR)/*

# For debugging variables
print-%:
	@echo $* = $($*)
