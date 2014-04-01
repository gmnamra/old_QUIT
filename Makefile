# C++11 settings are compiler specific
CXX_VER := $(shell cpp --version)
ifneq (,$(findstring clang,$(CXX_VER)))
	STDLIB      := -stdlib=libc++
	MOREFLAGS   := -Wno-deprecated-register
else
	THREADS     := -pthread
	STDLIB      := -lstdc++
	MOREFLAGS   :=
endif

# Pull in Eigen if necessary
EIGEN_DIR := eigen
EIGEN : $(EIGEN_DIR)/Eigen/Array
$(EIGEN_DIR)/Eigen/Array :
	@mkdir -p $(EIGEN_DIR)
	curl --location http://bitbucket.org/eigen/eigen/get/3.2.1.tar.gz > $(EIGEN_DIR).tar.gz
	tar --extract --file=$(EIGEN_DIR).tar.gz --strip-components=1 --directory=$(EIGEN_DIR)
	rm $(EIGEN_DIR).tar.gz

BUILD_DIR   := build
SOURCE_DIR  := src
INSTALL_DIR := ./bin

CXXFLAGS := -std=c++11 $(STDLIB) $(THREADS) -g -O3 -m64 -msse3 -mssse3 -msse4.1 -msse4.2 -Wfatal-errors -DAGILENT $(MOREFLAGS)
LDFLAGS  := -std=c++11 $(STDLIB) $(THREADS) -m64 -L$(BUILD_DIR)
INCLUDE    := -I$(EIGEN_DIR) -Isrc -Isrc/Agilent

#Rules for libNifti
NIFTI_DIR := Nifti
NIFTI_SRC := Nifti Internal ZipFile Extension
NIFTI_OBJ := $(patsubst %, $(BUILD_DIR)/$(NIFTI_DIR)/%.o, $(NIFTI_SRC))
$(BUILD_DIR)/$(NIFTI_DIR)/%.o : $(SOURCE_DIR)/$(NIFTI_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libNifti.a : $(NIFTI_OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $(NIFTI_OBJ)

#Rules for libAgilent
AGILENT_DIR := Agilent
AGILENT_SRC := fdf fdfFile FID FIDFile procpar util
AGILENT_OBJ := $(patsubst %, $(BUILD_DIR)/$(AGILENT_DIR)/%.o, $(AGILENT_SRC))
$(BUILD_DIR)/$(AGILENT_DIR)/%.o : $(SOURCE_DIR)/$(AGILENT_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libAgilent.a : $(AGILENT_OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $(AGILENT_OBJ)

#Rules for tools
TOOL_DIR   := Tools
TOOLS      := niihdr fdf2nii procparse
$(BUILD_DIR)/$(TOOL_DIR)/%.o : $(SOURCE_DIR)/$(TOOL_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(TOOLS)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(TOOL_DIR)/%.o | libAgilent.a libNifti.a
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz

#Rules for DESPOT
DESPOT      := afi despot1 despot1hifi despot2 despot2fm mcdespot mcsignal ssfpbands dixon phasemap
DESPOT_DIR  := DESPOT
DESPOT_SRC  := DESPOT DESPOT_Functors Model ThreadPool
DESPOT_TMPL := $(addprefix $(SOURCE_DIR)/$(DESPOT_DIR)/, RegionContraction.h DESPOT_Functors.h)
DESPOT_OBJ  := $(patsubst %, $(BUILD_DIR)/$(DESPOT_DIR)/%.o, $(DESPOT_SRC))
$(BUILD_DIR)/$(DESPOT_DIR)/%.o : $(SOURCE_DIR)/$(DESPOT_DIR)/%.cpp $(DESPOT_TMPL) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(DESPOT)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(DESPOT_DIR)/%.o $(DESPOT_OBJ) | libAgilent.a libNifti.a
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz

TARGETS := $(TOOLS) $(DESPOT)
LIB_TGT := libNifti.a libAgilent.a
all     : $(LIB_TGT) $(TARGETS)

clean :
	rm -rf $(BUILD_DIR)

install :
	@mkdir -p $(INSTALL_DIR)/
	cp $(addprefix $(BUILD_DIR)/, $(TARGETS)) $(INSTALL_DIR)/

.PHONY  : all install clean $(LIB_TGT) $(TARGETS) EIGEN
.DEFAULT_GOAL := all

$(TARGETS) : % : $(BUILD_DIR)/%
$(LIB_TGT) : % : $(BUILD_DIR)/%
