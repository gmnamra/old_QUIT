###################################################
#
# Makefile for DESPOT Fitting Programs
# Tobias Wood December 2012
#
###################################################

# Platform/system specific
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CPPPATH = /software/system/gcc/gcc-4.8.0
	LIBCPP = LD_RUN_PATH=$(CPPPATH)/lib64
	CXX = $(LIBCPP) $(CPPPATH)/bin/g++
	EIGEN  = /home/k1078535/Code/eigen
	INSTALL_DIR := /home/k1078535/Code
endif
ifeq ($(UNAME_S),Darwin)
	# Defaults work okay on Apple
	EIGEN  = /Users/Tobias/Code/eigen
	INSTALL_DIR := /Users/Tobias/Code/MR
endif

# Set up Paths
SRC_DIR = ../src
OBJ_DIR := obj
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib

# Set up all our compiler options
CXX_FLAGS = -g -std=c++11 -stdlib=libc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
INCLUDE = -I$(SRC_DIR) -I$(INSTALL_INC) -I$(EIGEN)
LD_FLAGS = -std=c++11 -stdlib=libc++ -m64 -O3 -L$(INSTALL_LIB)
LD_LIBS  = -lNifti -lz

# Top level build rules
.PHONY : all debug install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files

TARGETS = despot1 despot1hifi despot2 despot2fm mcdespot mcsignal phasemap afi

debug   : CXX_FLAGS+=-g
debug   : LD_FLAGS+=-g
debug   : all

all	: $(TARGETS)

install : $(INSTALL_BIN) 
	cp $(TARGETS) $(INSTALL_BIN)/

clean :
	rm -f $(TARGETS)
	rm -r $(OBJ_DIR)

%: $(OBJ_DIR)/%_main.o $(OBJ_DIR)/DESPOT.o $(OBJ_DIR)/DESPOT_Functors.o $(OBJ_DIR)/ThreadPool.o $(OBJ_DIR)/Model.o | $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) $(LD_LIBS)

$(OBJ_DIR)/mcdespot_main.o: $(SRC_DIR)/mcdespot_main.cpp $(SRC_DIR)/DESPOT_Functors.h $(SRC_DIR)/RegionContraction.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(OBJ_DIR)/mcsignal_main.o: $(SRC_DIR)/mcsignal_main.cpp $(SRC_DIR)/DESPOT_Functors.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<
	
$(OBJ_DIR)/despot2fm_main.o: $(SRC_DIR)/despot2fm_main.cpp $(SRC_DIR)/DESPOT_Functors.h $(SRC_DIR)/RegionContraction.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# For debugging variables
print-%:
	@echo $* = $($*)
