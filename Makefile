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
NIFTI_DIR  := Nifti
NIFTI_SRC  := Nifti Internal ZipFile Extension
NIFTI_OBJ  := $(addprefix $(BUILD_DIR)/$(NIFTI_DIR)/, $(addsuffix .o, $(NIFTI_SRC)))
NIFTI_HDR := $(addprefix $(SOURCE_DIR)/$(NIFTI_DIR)/, Nifti.h Nifti-inl.h ExtensionCodes.h)
$(BUILD_DIR)/$(NIFTI_DIR)/%.o : $(SOURCE_DIR)/$(NIFTI_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libNifti.a : $(NIFTI_OBJ) $(NIFTI_PHDR)
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

#Rules for libQUIT
QUIT_DIR   := QUIT
QUIT_SRC   := ThreadPool
QUIT_HDR   := $(addprefix $(SOURCE_DIR)/$(QUIT_DIR)/, Volume.h Volume-inl.h)
QUIT_OBJ   := $(addprefix $(BUILD_DIR)/$(QUIT_DIR)/, $(addsuffix .o, $(QUIT_SRC)))
$(BUILD_DIR)/$(QUIT_DIR)/%.o : $(SOURCE_DIR)/$(QUIT_DIR)/%.cpp | EIGEN libNifti.a
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libQUIT.a : $(QUIT_OBJ) $(QUIT_PHDR)
	@mkdir -p $(dir $@)
	ar rcs $@ $(QUIT_OBJ)

#Rules for tools
TOOL_DIR   := Tools
TOOLS      := niihdr procparse fdf2nii
PYTOOLS    := fdf2nii.py
$(BUILD_DIR)/$(TOOL_DIR)/%.o : $(SOURCE_DIR)/$(TOOL_DIR)/%.cpp $(NIFTI_PHDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(TOOLS)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(TOOL_DIR)/%.o | libAgilent.a libNifti.a
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) -lAgilent -lNifti -lz
$(addprefix $(BUILD_DIR)/, $(PYTOOLS)) :
	cp $(patsubst $(BUILD_DIR)/%, $(SOURCE_DIR)/$(TOOL_DIR)/%, $@) $(BUILD_DIR)

#Rules for DESPOT
DESPOT      := afi despot1 despot1hifi despot2 despot2fm mcdespot mcsignal ssfpbands phasemap
DESPOT_DIR  := DESPOT
DESPOT_SRC  := DESPOT DESPOT_Functors Model
DESPOT_HDR  := $(addprefix $(SOURCE_DIR)/$(DESPOT_DIR)/, RegionContraction.h DESPOT_Functors.h)
DESPOT_OBJ  := $(patsubst %, $(BUILD_DIR)/$(DESPOT_DIR)/%.o, $(DESPOT_SRC))
$(BUILD_DIR)/$(DESPOT_DIR)/%.o : $(SOURCE_DIR)/$(DESPOT_DIR)/%.cpp $(DESPOT_HDR) $(QUIT_HDR) $(NIFTI_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(DESPOT)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(DESPOT_DIR)/%.o $(DESPOT_OBJ) | libAgilent.a libNifti.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz

#Rules for Misc
MISC     := dixon
MISC_DIR := Misc
$(BUILD_DIR)/$(MISC_DIR)/%.o : $(SOURCE_DIR)/$(MISC_DIR)/%.cpp $(QUIT_HDR) $(NIFTI_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(MISC)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(MISC_DIR)/%.o $(MISC_OBJ) | libAgilent.a libNifti.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz

TARGETS := $(TOOLS) $(PYTOOLS) $(DESPOT) $(MISC)
LIB_TGT := libNifti.a libAgilent.a libQUIT.a
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
