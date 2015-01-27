# C++11 settings are compiler specific
CXX_VER := $(shell cpp --version)
ifneq (,$(findstring clang,$(CXX_VER)))
	STDLIB      := -stdlib=libc++
	MOREFLAGS   := -Wno-deprecated-register
else
	THREADS     := -pthread
	STDLIB      := -lstdc++
	MOREFLAGS   := -Wno-narrowing
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

$(BUILD_DIR) :
	@mkdir -p $(BUILD_DIR)
$(INSTALL_DIR) :
	@mkdir -p $(INSTALL_DIR)

CXXFLAGS := -std=c++11 $(STDLIB) $(THREADS) -march=native -Wfatal-errors $(MOREFLAGS)
LDFLAGS  := -std=c++11 $(STDLIB) $(THREADS) -march=native -L$(BUILD_DIR)
INCLUDE    := -I$(EIGEN_DIR) -Isrc -Isrc/Agilent

#Rules for libNifti
NIFTI_DIR  := Nifti
NIFTI_SRC  := Nifti Header Internal ZipFile Extension
NIFTI_OBJ  := $(addprefix $(BUILD_DIR)/$(NIFTI_DIR)/, $(addsuffix .o, $(NIFTI_SRC)))
NIFTI_HDR  := $(addprefix $(SOURCE_DIR)/$(NIFTI_DIR)/, Nifti.h Nifti-inl.h Enum.h ExtensionCodes.h)
$(BUILD_DIR)/$(NIFTI_DIR)/%.o : $(SOURCE_DIR)/$(NIFTI_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libNifti.a : $(NIFTI_OBJ) $(NIFTI_HDR)
	@mkdir -p $(dir $@)
	ar rcs $@ $(NIFTI_OBJ)

#Rules for libAgilent
AGILENT_DIR := Agilent
AGILENT_SRC := fdf fdfFile FID FIDFile procpar util
AGILENT_OBJ := $(patsubst %, $(BUILD_DIR)/$(AGILENT_DIR)/%.o, $(AGILENT_SRC))
AGILENT_HDR := $(addprefix $(SOURCE_DIR)/$(AGILENT_DIR)/, $(addsuffix .h, $(AGILENT_SRC)))
$(BUILD_DIR)/$(AGILENT_DIR)/%.o : $(SOURCE_DIR)/$(AGILENT_DIR)/%.cpp | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libAgilent.a : $(AGILENT_OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $(AGILENT_OBJ)

#Rules for libQUIT
QUIT_DIR   := QUIT
QUIT_SRC   := Util ThreadPool
QUIT_TPL   := MultiArray MultiArray-inl Volume Volume-inl
QUIT_HDR   := $(addprefix $(SOURCE_DIR)/$(QUIT_DIR)/, $(addsuffix .h, QUIT $(QUIT_SRC) $(QUIT_TPL)))
QUIT_OBJ   := $(addprefix $(BUILD_DIR)/$(QUIT_DIR)/, $(addsuffix .o, $(QUIT_SRC)))
$(BUILD_DIR)/$(QUIT_DIR)/%.o : $(SOURCE_DIR)/$(QUIT_DIR)/%.cpp $(NIFTI_HDR) $(AGILENT_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libQUIT.a : $(QUIT_OBJ) $(QUIT_HDR)
	@mkdir -p $(dir $@)
	ar rcs $@ $(QUIT_OBJ)

#Rules for tools
TOOL_DIR   := Tools
TOOLS      := niihdr niiext niicreate niicomplex niinudge niigrad procparse fdf2nii fdfbval calctfm
PYTOOLS    := fdf2nii.py
$(BUILD_DIR)/$(TOOL_DIR)/%.o : $(SOURCE_DIR)/$(TOOL_DIR)/%.cpp $(NIFTI_HDR) $(QUIT_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(TOOLS)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(TOOL_DIR)/%.o libNifti.a libAgilent.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $< -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz
$(addprefix $(BUILD_DIR)/, $(PYTOOLS)) : $(BUILD_DIR} | $(TOOL_DIR)/$@
	cp $(patsubst $(BUILD_DIR)/%, $(SOURCE_DIR)/$(TOOL_DIR)/%, $@) $(BUILD_DIR)/

#Rules for DESPOT
DESPOT      := afi despot1 despot1hifi despot2 despot2fm mcdespot mcsignal ssfpbands phasemap t2star
DESPOT_DIR  := DESPOT
DESPOT_SRC  := DESPOT DESPOT_Functors Model
DESPOT_HDR  := $(addprefix $(SOURCE_DIR)/$(DESPOT_DIR)/, RegionContraction.h DESPOT_Functors.h DESPOT.h Model.h)
DESPOT_OBJ  := $(patsubst %, $(BUILD_DIR)/$(DESPOT_DIR)/%.o, $(DESPOT_SRC))
$(BUILD_DIR)/$(DESPOT_DIR)/%.o : $(SOURCE_DIR)/$(DESPOT_DIR)/%.cpp $(DESPOT_HDR) $(QUIT_HDR) $(NIFTI_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(DESPOT)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(DESPOT_DIR)/%.o $(DESPOT_OBJ) libNifti.a libAgilent.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $< $(DESPOT_OBJ) -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz

#Rules for Misc
MISC     := dixon
MISC_DIR := Misc
$(BUILD_DIR)/$(MISC_DIR)/%.o : $(SOURCE_DIR)/$(MISC_DIR)/%.cpp $(QUIT_HDR) $(NIFTI_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(MISC)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(MISC_DIR)/%.o $(MISC_OBJ) libNifti.a libAgilent.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $< -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz

TARGETS := $(TOOLS) $(PYTOOLS) $(DESPOT) $(MISC)
LIB_TGT := libNifti.a libAgilent.a libQUIT.a

all     : CXXFLAGS += -O3
all     : $(LIB_TGT) $(TARGETS)

debug   : CXXFLAGS += -g
debug   : LDFLAGS  += -g
debug   : $(LIB_TGT) $(TARGETS)

clean :
	rm -rf $(BUILD_DIR)

install :
	@mkdir -p $(INSTALL_DIR)/
	cp $(addprefix $(BUILD_DIR)/, $(TARGETS)) $(INSTALL_DIR)/
	chmod ugo+rx $(addprefix $(INSTALL_DIR)/, $(TARGETS))

.PHONY  : all install clean
.INTERMEDIATE : $(LIB_TGT) $(TARGETS) EIGEN
.DEFAULT_GOAL := all

$(TARGETS) : % : $(BUILD_DIR)/%
$(LIB_TGT) : % : $(BUILD_DIR)/%
