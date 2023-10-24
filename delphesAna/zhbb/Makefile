#ROOT_DIR = /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.14.04-0d8dc/x86_64-slc6-gcc7-opt/
ROOT_DIR = `root-config --prefix`
# Delphes directory                                                                                                                                                                                                                                                                                                                                                         
#DELPHES_DIR = ../..
DELPHES_DIR=${DELPHES_HOME}

CXXFLAGS += $(shell $(ROOT_DIR)/bin/root-config --cflags)
CXXFLAGS += -I$(DELPHES_DIR) -I$(DELPHES_DIR)/external

LDFLAGS += $(shell $(ROOT_DIR)/bin/root-config --libs) -lEG
LDFLAGS += -Wl,-rpath,$(DELPHES_DIR) -L$(DELPHES_DIR) -lDelphes

all: zhbb_analyze
clean:
	rm -f zhbb_analyze

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
