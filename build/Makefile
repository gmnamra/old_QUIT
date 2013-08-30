###################################################
#
# Makefile for libAgilent
# Tobias Wood December 2012
#
###################################################

# Set up variables
SRC_DIR = ../src
OBJ_DIR := obj
INSTALL_DIR := /home/k1078535/Code
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/system/gcc/gcc-4.8.0/lib64
EIGEN = /home/k1078535/Code/eigen
TARGETS = libAgilent.a procparse fdf2nii

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files
all : $(TARGETS)

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/g++
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -pthread -Wfatal-errors $(DEBUG)
AR = LD_RUN_PATH=$(LIBCPP) ar rcs
LD_FLAGS = -std=c++11 -lstdc++ -O3 -pthread -L. -L$(INSTALL_LIB)
INCLUDE = -I$(SRC_DIR) -I$(EIGEN) -I$(INSTALL_INC)

# Source files for libAgilent
AG_BASE = util procpar fdf fdfFile FID FIDFile
AG_HDR = $(patsubst %, $(SRC_DIR)/%.h, $(AG_BASE))
AG_OBJ = $(patsubst %, $(OBJ_DIR)/%.o, $(AG_BASE))

# Build rules
%: $(OBJ_DIR)/%_main.o libAgilent.a | $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) -lAgilent -lNifti -lz

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

libAgilent.a : $(AG_OBJ) | $(OBJ_DIR)
	$(AR) libAgilent.a $(AG_OBJ)

install : | $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp procparse fdf2nii $(INSTALL_BIN)/
	cp libAgilent.a $(INSTALL_LIB)/
	cp $(AG_HDR) $(INSTALL_INC)/

clean :
	rm -f $(TARGETS)
	rm -r $(OBJ_DIR)

# For debugging variables
print-%:
	@echo $* = $($*)
