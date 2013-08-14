###################################################
#
# Makefile for DESPOT Fitting Programs
# Tobias Wood December 2012
#
###################################################

# Set up Paths
SRC_DIR = ..
OBJ_DIR := obj
INSTALL_DIR := /home/k1078535/Code
INSTALL_BIN = $(INSTALL_DIR)/bin
INSTALL_INC = $(INSTALL_DIR)/include
INSTALL_LIB = $(INSTALL_DIR)/lib
LIBCPP = /software/system/gcc/gcc-4.8.0/lib64
EIGEN  = /home/k1078535/Code/eigen

# Set up all our compiler options
CXX = LD_RUN_PATH=$(LIBCPP) /software/system/gcc/gcc-4.8.0/bin/g++
CXX_FLAGS = -std=c++11 -lstdc++ -m64 -msse3 -mssse3 -msse4.1 -msse4.2 -pthread -Wfatal-errors $(DEBUG) -DHAVE_NRECON
INCLUDE = -I$(SRC_DIR) -I$(INSTALL_INC) -I$(EIGEN)
LD_FLAGS = -std=c++11 -lstdc++ -O3 -pthread -L$(INSTALL_LIB)
LD_LIBS  = -lrecon -lNifti -lz
DEBUG_FLAGS = -g
RELEASE_FLAGS = -O3

# Top level build rules
.PHONY : release debug install clean
.PRECIOUS: $(OBJ_DIR)/%.o # Stop make removing the object files

TARGETS = despot1 despot-hifi despot2 despot2-fm mcdespot mcfinite phasemap afi

release : CXX_FLAGS+=$(RELEASE_FLAGS)
release : $(TARGETS)

debug   : CXX_FLAGS+=$(DEBUG_FLAGS)
debug   : $(TARGETS)

install : $(INSTALL_BIN) 
	cp $(TARGETS) $(INSTALL_BIN)/

clean :
	rm -f $(TARGETS)
	rm -r $(OBJ_DIR)

%: $(OBJ_DIR)/%_main.o $(OBJ_DIR)/DESPOT.o $(OBJ_DIR)/ThreadPool.o | $(OBJ_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) $(LD_LIBS)

$(OBJ_DIR)/mcfinite_main.o: $(SRC_DIR)/mcdespot_main.cpp $(SRC_DIR)/DESPOT_Functors.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -DUSE_MCFINITE -o $@ $<

$(OBJ_DIR)/mcdespot_main.o: $(SRC_DIR)/mcdespot_main.cpp $(SRC_DIR)/DESPOT_Functors.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(OBJ_DIR)/despot2-fm_main.o: $(SRC_DIR)/despot2-fm_main.cpp $(SRC_DIR)/DESPOT_Functors.h | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

# Create directories
$(OBJ_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# For debugging variables
print-%:
	@echo $* = $($*)
