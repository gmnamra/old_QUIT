###################################################
#
# Makefile for libNifti
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
SRC_DIR     := src/Nifti
INC_DIR     := src
BLD_DIR     := build
INSTALL_BIN := $(INSTALL_DIR)/bin
INSTALL_INC := $(INSTALL_DIR)/include
INSTALL_LIB := $(INSTALL_DIR)/lib

# Set up all our compiler options
CXX_FLAGS := -std=c++11 $(STDLIB) $(THREADS) -m64 -g -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors -DAGILENT $(DEBUG)
INCLUDE   := -I$(INC_DIR) -I$(INSTALL_INC) -I$(EIGEN)
LD_FLAGS  := -std=c++11 $(STDLIB) $(THREADS) -m64 -g -L$(INSTALL_LIB)
LD_LIBS   := -lAgilent -lNifti -lz

# Create directories
$(BLD_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Build rules for libNifti
NIFTI_BASE = Nifti Internal ZipFile Extension
NIFTI_OBJ  = $(patsubst %, $(BLD_DIR)/%.o, $(NIFTI_BASE))
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<
libNifti.a : $(NIFTI_OBJ) | $(BLD_DIR)
	ar rcs libNifti.a $(NIFTI_OBJ)

# Build rules for utils (niihdr only at the moment)
$(BLD_DIR)/niihdr.o: src/niihdr.cpp | libNifti.a $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

niihdr : $(BLD_DIR)/niihdr.o | libNifti.a $(BLD_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) -lNifti -lz

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(BLD_DIR)/%.o # Stop make removing the object files

TARGETS = libNifti.a niihdr

all     : $(TARGETS)

install : | $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	mkdir -p $(INSTALL_INC)/Nifti
	cp $(SRC_DIR)/Nifti.h $(SRC_DIR)/ExtensionCodes.h $(SRC_DIR)/Nifti-inl.h $(SRC_DIR)/ZipFile.h $(SRC_DIR)/Volume.h $(SRC_DIR)/Volume-inl.h $(INSTALL_INC)/Nifti/
	cp libNifti.a $(INSTALL_LIB)/
	cp niihdr $(INSTALL_BIN)/

clean : 
	rm -r $(BLD_DIR)
	rm -f $(TARGETS)

# For debugging variables
print-%:
	@echo $* = $($*)
