v###################################################
#
# Makefile for libNifti
# Tobias Wood December 2012
#
###################################################

# Set up Paths
SRC_DIR = ../src
OBJ_DIR := obj
INSTALL_DIR := .
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/system/gcc/gcc-4.8.0/lib
EIGEN  = /software/local/k1078535/include/eigen3

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/gcc
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
AR = LD_RUN_PATH=$(LIBCPP) ar rcs
LD_FLAGS = -O3 -L$(OBJ_DIR)
INCLUDE = -I$(SRC_DIR) -I$(EIGEN)

# Source files for libNifti

NIFTI_BASE = Nifti
NIFTI_OBJ  = $(patsubst %, $(OBJ_DIR)/%.o, $(NIFTI_BASE))
NIFTI_HDR  = $(patsubst %, $(SRC_DIR)/%.h, $(NIFTI_BASE)) $(SRC_DIR)/nifti1.h

# General Build Rule
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

libNifti : $(NIFTI_OBJ) | $(OBJ_DIR)
	$(AR) $(LIB_OPTIONS) $(OBJ_DIR)/libNifti.a $(NIFTI_OBJ)

nifti_hdr : $(OBJ_DIR)/nifti_hdr.o libNifti | $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) -LNifti

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files

TARGETS = libNifti nifti_hdr

all     : $(TARGETS)

install : | $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp $(NIFTI_HDR) $(INSTALL_INC)/
	cp $(OBJ_DIR)/libNifti.a $(INSTALL_LIB)/
	cp $(OBJ_DIR)/nifti_hdr $(INSTALL_BIN)/

clean : 
	rm -r $(OBJ_DIR)
	rm -f $(TARGETS)

# For debugging variables
print-%:
	@echo $* = $($*)
