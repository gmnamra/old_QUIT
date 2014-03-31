# C++11 settings are compiler specific
CXX_VER := $(shell cpp --version)
ifneq (,$(findstring clang,$(CXX_VER)))
	STDLIB      := -stdlib=libc++
	EIGEN       := $(HOME)/Code/eigen
	INSTALL_DIR := $(HOME)/Code/MR
	MOREFLAGS   := -Wno-deprecated-register
else
	THREADS     := -pthread
	STDLIB      := -lstdc++
	INSTALL_DIR := $(HOME)/Code
	EIGEN       := $(INSTALL_DIR)/eigen
	MOREFLAGS   :=
endif

BUILD_DIR  := build
SOURCE_DIR := src
INCLUDE    := -I$(EIGEN) -Isrc -Isrc/Agilent

CXXFLAGS := -std=c++11 $(STDLIB) $(THREADS) -g -O3 -m64 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors -DAGILENT $(MOREFLAGS)
LDFLAGS  := -std=c++11 $(STDLIB) $(THREADS) -m64 -L$(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $@

DESPOT  := afi despot1 despot1hifi despot2 despot2fm mcdespot mcsignal ssfpbands dixon
TARGETS := procparse fdf2nii $(DESPOT)
LIB_TGT := libNifti.a libAgilent.a
all     : $(LIB_TGT) $(TARGETS)
.PHONY  : all $(LIB_TGT) $(TARGETS)
.DEFAULT_GOAL := all

$(TARGETS) : % : $(BUILD_DIR)/%
$(LIB_TGT) : % : $(BUILD_DIR)/%

#Rules for libNifti
NIFTI_DIR := Nifti
NIFTI_SRC := Nifti Internal ZipFile Extension
NIFTI_OBJ := $(patsubst %, $(BUILD_DIR)/$(NIFTI_DIR)/%.o, $(NIFTI_SRC))
$(BUILD_DIR)/$(NIFTI_DIR) :
	mkdir -p $@
$(BUILD_DIR)/$(NIFTI_DIR)/%.o : $(SOURCE_DIR)/$(NIFTI_DIR)/%.cpp | $(BUILD_DIR)/$(NIFTI_DIR)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libNifti.a : $(NIFTI_OBJ) | $(BUILD_DIR)/$(NIFTI_DIR)
	ar rcs $@ $(NIFTI_OBJ)

#Rules for libAgilent
AGILENT_DIR := Agilent
$(BUILD_DIR)/$(AGILENT_DIR):
	mkdir -p $@
AGILENT_SRC := fdf fdfFile FID FIDFile procpar util
AGILENT_OBJ := $(patsubst %, $(BUILD_DIR)/$(AGILENT_DIR)/%.o, $(AGILENT_SRC))
$(BUILD_DIR)/$(AGILENT_DIR) :
	mkdir -p $@
$(BUILD_DIR)/$(AGILENT_DIR)/%.o : $(SOURCE_DIR)/$(AGILENT_DIR)/%.cpp | $(BUILD_DIR)/$(AGILENT_DIR)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libAgilent.a : $(AGILENT_OBJ) | $(BUILD_DIR)
	ar rcs $@ $(AGILENT_OBJ)

#Rules for tools
TOOL_DIR := Tools
$(BUILD_DIR)/$(TOOL_DIR):
	mkdir -p $@
#TOOL_SRC := procparse_main fdf2nii_main niihdr
#TOOL_OBJ := $(patsubst %, $(BUILD_DIR)/$(TOOL_DIR)/%.o, $(TOOL_SRC))
$(BUILD_DIR)/$(TOOL_DIR)/%.o : $(SOURCE_DIR)/$(TOOL_DIR)/%.cpp | $(BUILD_DIR)/$(TOOL_DIR)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/procparse : $(BUILD_DIR)/$(TOOL_DIR)/procparse.o | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/fdf2nii : $(BUILD_DIR)/$(TOOL_DIR)/fdf2nii.o | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/niihdr : $(BUILD_DIR)/$(TOOL_DIR)/niihdr.o | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lNifti -lz

#Rules for DESPOT
DESPOT_DIR := DESPOT
$(BUILD_DIR)/$(DESPOT_DIR):
	mkdir -p $@
DESPOT_SRC := DESPOT DESPOT_Functors Model ThreadPool
DESPOT_TEMPLATES := RegionContraction.h DESPOT_Functors.h
DESPOT_OBJ := $(patsubst %, $(BUILD_DIR)/$(DESPOT_DIR)/%.o, $(DESPOT_SRC))
$(BUILD_DIR)/$(DESPOT_DIR)/%.o : $(SOURCE_DIR)/$(DESPOT_DIR)/%.cpp | $(BUILD_DIR)/$(DESPOT_DIR)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/afi : $(BUILD_DIR)/$(DESPOT_DIR)/afi_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/despot1 : $(BUILD_DIR)/$(DESPOT_DIR)/despot1_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/despot1hifi : $(BUILD_DIR)/$(DESPOT_DIR)/despot1hifi_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/despot2 : $(BUILD_DIR)/$(DESPOT_DIR)/despot2_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/despot2fm : $(BUILD_DIR)/$(DESPOT_DIR)/despot2fm_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/mcdespot : $(BUILD_DIR)/$(DESPOT_DIR)/mcdespot_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/mcsignal : $(BUILD_DIR)/$(DESPOT_DIR)/mcsignal_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/ssfpbands : $(BUILD_DIR)/$(DESPOT_DIR)/ssfpbands_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/dixon : $(BUILD_DIR)/$(DESPOT_DIR)/dixon_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(BUILD_DIR)/phasemap : $(BUILD_DIR)/$(DESPOT_DIR)/phasemap_main.o $(DESPOT_OBJ) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
