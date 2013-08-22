###################################################
#
# Makefile for librecon
# Tobias Wood December 2012
#
###################################################

# Set up Paths
SRC_DIR = ../src
BUILD_DIR := obj
INSTALL_DIR := .
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/system/gcc/gcc-4.8.0/lib64
EIGEN = ~/Code/eigen

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/g++
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -pthread -Wfatal-errors $(DEBUG)
INCLUDE = -I$(SRC_DIR) -I$(EIGEN)
LD_FLAGS = -std=c++11 -lstdc++ -O3 -pthread -L$(INSTALL_LIB)
LD_LIBS  = 

# Source files for libagilent
AG_BASE = FID FIDFile recon_util
AG_HDR = $(patsubst %, $(SRC_DIR)/%.h, $(AG_BASE))
AG_OBJ = $(patsubst %, $(OBJ_DIR)/%.o, $(AG_BASE))

# Build rules
%: $(OBJ_DIR)/%_main.o | $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) $(LD_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

libagilent : $(AG_OBJ) | $(OBJ_DIR)
	$(AR) libagilent.a $(RECON_DEPS)

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files

TARGETS = libagilent procparse

all	: $(TARGETS)

install : $(INSTALL_BIN) 
	cp procparse $(INSTALL_BIN)/
	cp libagilent $(INSTALL_LIB)/
	cp $(AG_HDR) $(INSTALL_HDR)/

clean :
	rm -f $(TARGETS)
	rm -r $(OBJ_DIR)

# For debugging variables
print-%:
	@echo $* = $($*)
