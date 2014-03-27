###################################################
#
# Makefile for DESPOT Fitting Programs
# Tobias Wood December 2012
#
###################################################

# C++11 settings are compiler specific
CXX_VER := $(shell cpp --version)
ifneq (,$(findstring clang,$(CXX_VER)))
	STDLIB      := -stdlib=libc++
	EIGEN       := $(HOME)/Code/eigen
	INSTALL_DIR := $(HOME)/Code/MR
else
	THREADS     := -pthread
	STDLIB      := -lstdc++
	INSTALL_DIR := $(HOME)/Code
	EIGEN       := $(INSTALL_DIR)/eigen
endif

# Set up Paths
SRC_DIR     := src
BLD_DIR     := build
INSTALL_BIN := $(INSTALL_DIR)/bin
INSTALL_INC := $(INSTALL_DIR)/include
INSTALL_LIB := $(INSTALL_DIR)/lib

# Set up all our compiler options
CXX_FLAGS := -std=c++11 $(STDLIB) $(THREADS) -m64 -g -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors -DAGILENT $(DEBUG)
INCLUDE   := -I$(SRC_DIR) -I$(INSTALL_INC) -I$(EIGEN)
LD_FLAGS  := -std=c++11 $(STDLIB) $(THREADS) -m64 -g -L$(INSTALL_LIB)
LD_LIBS   := -lAgilent -lNifti -lz

# Top level build rules
.PHONY : all debug install clean
.PRECIOUS: $(BLD_DIR)/%.o # Stop make removing the object files

%: $(BLD_DIR)/%_main.o $(BLD_DIR)/DESPOT.o $(BLD_DIR)/DESPOT_Functors.o $(BLD_DIR)/ThreadPool.o $(BLD_DIR)/Model.o | $(BLD_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) $(LD_LIBS)

$(BLD_DIR)/mcdespot_main.o: $(SRC_DIR)/mcdespot_main.cpp $(SRC_DIR)/DESPOT_Functors.h $(SRC_DIR)/RegionContraction.h | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(BLD_DIR)/mcsignal_main.o: $(SRC_DIR)/mcsignal_main.cpp $(SRC_DIR)/DESPOT_Functors.h | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<
	
$(BLD_DIR)/despot2fm_main.o: $(SRC_DIR)/despot2fm_main.cpp $(SRC_DIR)/DESPOT_Functors.h $(SRC_DIR)/RegionContraction.h | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

# Create directories
$(BLD_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

TARGETS = despot1 despot1hifi despot2 despot2fm mcdespot mcsignal phasemap afi ssfpbands dixon

all	: $(TARGETS)

install : $(INSTALL_BIN)
	cp $(TARGETS) $(INSTALL_BIN)/

release : CXX_FLAGS += -O3
release : LD_FLAGS += -O3
release : all
clean :
	rm -f $(TARGETS)
	rm -f $(BLD_DIR)/*.o
# For debugging variables
print-%:
	@echo $* = $($*)
