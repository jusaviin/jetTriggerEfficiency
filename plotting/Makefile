PROGRAM       = DrawingClasses
version       = EEC
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -D$(version) # Add -g for debugger
LD            = g++
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

CXXFLAGS  += $(INCS)

HDRS	+= TriggerHistogramManager.h TriggerCard.h AlgorithmLibrary.h
DEPENDENCIES += ../src/TriggerHistograms.o ../src/ConfigurationCard.o

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM) :    $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -shared -lEG -L$(PWD) $(CXXFLAGS) $(OBJS) $(DEPENDENCIES) $(LDFLAGS) -o $(PROGRAM).so
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) $(PROGRAM).o $(PROGRAM).so

cl:  clean $(PROGRAM)
