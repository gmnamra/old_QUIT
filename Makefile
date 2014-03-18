###################################################
#
# Makefile for libAgilent
# Tobias Wood December 2012
#
###################################################

# Platform/system specific
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

# Set up variables
SRC_DIR := src
PYTHON_DIR := python
BLD_DIR := build
INSTALL_BIN := $(INSTALL_DIR)/bin
INSTALL_INC := $(INSTALL_DIR)/include
INSTALL_LIB := $(INSTALL_DIR)/lib

TARGETS = libAgilent.a procparse fdf2nii fdfhdr

# Top level build rules
.PHONY : all install clean
.PRECIOUS: $(BLD_DIR)/%.o # Stop make removing the object files
all : $(TARGETS)

# Create directories
$(BLD_DIR) $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB):
	mkdir -p $@

# Set up all our compiler options
CXX_FLAGS := -g -std=c++11 -stdlib=libc++ -m64 -O3 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors $(DEBUG)
LD_FLAGS  := -g -std=c++11 -stdlib=libc++ -m64 -O3 -L. -L$(INSTALL_LIB)
INCLUDE   := -I$(SRC_DIR) -I$(EIGEN) -I$(INSTALL_INC)

# Source files for libAgilent
AG_BASE := util procpar fdf fdfFile FID FIDFile
AG_HDR  := $(patsubst %, $(SRC_DIR)/%.h, $(AG_BASE))
AG_OBJ  := $(patsubst %, $(BLD_DIR)/%.o, $(AG_BASE))

# Build rules
%: $(BLD_DIR)/%_main.o libAgilent.a | $(BLD_DIR)
	$(CXX) $^ -o $@ $(LD_FLAGS) -lAgilent -lNifti -lz

$(BLD_DIR)/%.o : $(SRC_DIR)/%.cpp | $(BLD_DIR)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<

libAgilent.a : $(AG_OBJ) | $(BLD_DIR)
	$(LIBCPP) ar rcs libAgilent.a $(AG_OBJ)

install : | $(INSTALL_BIN) $(INSTALL_INC) $(INSTALL_LIB)
	cp procparse fdf2nii fdfhdr $(PYTHON_DIR)/fdf2nii.py $(INSTALL_BIN)/
	cp libAgilent.a $(INSTALL_LIB)/
	cp $(AG_HDR) $(INSTALL_INC)/

clean :
	rm -f $(TARGETS)
	rm -r $(BLD_DIR)

# For debugging variables
print-%:
	@echo $* = $($*)
