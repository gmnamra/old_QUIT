###################################################
#
# Makefile for libNifti
# Tobias Wood December 2012
#
###################################################

# Platform/system specific
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CPPPATH     := /software/system/gcc/gcc-4.8.0
	LIBCPP      := LD_RUN_PATH=$(CPPPATH)/lib64
	CXX         := $(LIBCPP) $(CPPPATH)/bin/g++
	EIGEN       := /home/k1078535/Code/eigen
	INSTALL_DIR := /home/k1078535/Code
endif
ifeq ($(UNAME_S),Darwin)
	# Defaults work okay on Apple
	EIGEN       := /Users/Tobias/Code/eigen
	INSTALL_DIR := /Users/Tobias/Code/MR
endif

# Set up Paths
SRC_DIR := src/Nifti
INC_DIR := src
OBJ_DIR := build
INSTALL_DIR := /Users/Tobias/Code/MR
INSTALL_BIN := $(INSTALL_DIR)/bin
INSTALL_INC := $(INSTALL_DIR)/include
INSTALL_LIB := $(INSTALL_DIR)/lib

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Set up all our compiler options
CXX_FLAGS = -g -std=c++11 -stdlib=libc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
LD_FLAGS = -g -std=c++11 -stdlib=libc++ -m64 -O3 -L.
INCLUDE = -I$(INC_DIR) -I$(EIGEN)

# Build rules for libNifti
NIFTI_BASE = Nifti Internal ZipFile Extension
NIFTI_OBJ  = $(patsubst %, $(OBJ_DIR)/%.o, $(NIFTI_BASE))
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<
libNifti.a : $(NIFTI_OBJ) | $(OBJ_DIR)
	$(LIBCPP) ar rcs libNifti.a $(NIFTI_OBJ)

# Build rules for utils (nifti_hdr only at the moment)
$(OBJ_DIR)/nifti_hdr.o: src/nifti_hdr.cpp | libNifti.a $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

nifti_hdr : $(OBJ_DIR)/nifti_hdr.o | libNifti.a $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) -lNifti -lz

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files

TARGETS = libNifti.a nifti_hdr

all     : $(TARGETS)

install : | $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	mkdir -p $(INSTALL_INC)/Nifti
	cp $(SRC_DIR)/Nifti.h $(SRC_DIR)/ExtensionCodes.h $(SRC_DIR)/Nifti-inl.h $(SRC_DIR)/ZipFile.h $(INSTALL_INC)/Nifti/
	cp libNifti.a $(INSTALL_LIB)/
	cp nifti_hdr $(INSTALL_BIN)/

clean : 
	rm -r $(OBJ_DIR)
	rm -f $(TARGETS)

# For debugging variables
print-%:
	@echo $* = $($*)
