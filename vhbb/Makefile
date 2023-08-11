# 'make' builds everything
# 'make clean' deletes everything except source files and Makefile
#
# You need to set ROOT_DIR, DELPHES_DIR and FASTJET_DIR for your project.
#
# ROOT directory
#ROOT_DIR = /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc7-opt/
ROOT_DIR = /Users/gaetano/programs/root_file/root_build_ventura_root_v6.26_00_arm_xcode14_2/
# Delphes directory
DELPHES_DIR = /Users/gaetano/Documents/universita/SnowMass2020/Analysis/spandansgeneration/MG5_aMC_v3_4_2/Delphes

#

CXXFLAGS += $(shell $(ROOT_DIR)/bin/root-config --cflags)
CXXFLAGS += -I$(DELPHES_DIR) -I$(DELPHES_DIR)/external

LDFLAGS += $(shell $(ROOT_DIR)/bin/root-config --libs) -lEG
LDFLAGS += -Wl,-rpath,$(DELPHES_DIR) -L$(DELPHES_DIR) -lDelphes

all: zhbb_analyze

clean:
	rm -f zhbb_analyze

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
