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
	curl --location http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz > $(EIGEN_DIR).tar.gz
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

#Rules for QUIT & libQUIT
QUIT_DIR   := QUIT
LIB_SRC   := Util ThreadPool
LIB_TPL   := MultiArray MultiArray-inl Volume Volume-inl
LIB_HDR   := $(addprefix $(SOURCE_DIR)/$(QUIT_DIR)/, $(addsuffix .h, QUIT $(LIB_SRC) $(LIB_TPL)))
LIB_OBJ   := $(addprefix $(BUILD_DIR)/$(QUIT_DIR)/, $(addsuffix .o, $(LIB_SRC)))
EXE      := afi despot1 despot1hifi despot2 despot2fm mcdespot mcsignal ssfpbands phasemap multiecho dixon
EXE_SRC  := SignalEquations Models Sequence
EXE_HDR  := $(addprefix $(SOURCE_DIR)/$(QUIT_DIR)/, SignalEquations.h Models.h Sequence.h RegionContraction.h)
EXE_OBJ  := $(patsubst %, $(BUILD_DIR)/$(QUIT_DIR)/%.o, $(EXE_SRC))
$(BUILD_DIR)/$(QUIT_DIR)/%.o : $(SOURCE_DIR)/$(QUIT_DIR)/%.cpp $(QUIT_HDR) $(LIB_HDR) $(NIFTI_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(BUILD_DIR)/libQUIT.a : $(LIB_OBJ) $(LIB_HDR)
	@mkdir -p $(dir $@)
	ar rcs $@ $(LIB_OBJ)
$(addprefix $(BUILD_DIR)/, $(EXE)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(QUIT_DIR)/%.o $(EXE_OBJ) libNifti.a libAgilent.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $< $(EXE_OBJ) -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz

#Rules for Tools
TOOL_DIR   := Tools
TOOLS      := niihdr niiext niicreate niicomplex niinudge niigrad procparse fdf2nii fdfbval calctfm
PYTOOLS    := fdf2nii.py
$(BUILD_DIR)/$(TOOL_DIR)/%.o : $(SOURCE_DIR)/$(TOOL_DIR)/%.cpp $(NIFTI_HDR) $(LIB_HDR) | EIGEN
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<
$(addprefix $(BUILD_DIR)/, $(TOOLS)) : $(BUILD_DIR)/% : $(BUILD_DIR)/$(TOOL_DIR)/%.o libNifti.a libAgilent.a libQUIT.a
	@mkdir -p $(dir $@)
	$(CXX) $< -o $@ $(LDFLAGS) -lQUIT -lAgilent -lNifti -lz
$(addprefix $(BUILD_DIR)/, $(PYTOOLS)) : $(BUILD_DIR)
	cp $(patsubst $(BUILD_DIR)/%, $(SOURCE_DIR)/$(TOOL_DIR)/%, $@) $(BUILD_DIR)/

TARGETS := $(EXE) $(TOOLS) $(PYTOOLS)
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
