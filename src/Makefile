CXX=g++
CXXFLAGS=-Wall -fopenmp
IDIR=

OBJDIR_RELEASE = ../build/release
OBJDIR_DEBUG = ../build/debug
LIBDIR_RELEASE = ../lib/release
LIBDIR_DEBUG = ../lib/debug

debug: CXXFLAGS += -DDEBUG -g -O0
debug: OBJDIR = $(OBJDIR_DEBUG)
debug: LIBDIR = $(LIBDIR_DEBUG)
debug: DIR=debug
debug: ../bin/debug/match

release: CXXFLAGS += -O3 -DNDEBUG
debug: OBJDIR = $(OBJDIR_RELEASE)
debug: LIBDIR = $(LIBDIR_RELEASE)
release: DIR=release
release: ../bin/release/match

testDatabase: ../test/bin/testDatabase

srcs=match.cc\
     globals.cc \
     dp.cc \
     Scorer.cc \
	 utils.cc \
	 debugUtils.cc \
	 ContigMapData.cc \
	 MatchResult.cc \
	 OpticalMapData.cc \
	 pugixml.cc \
	 StandardMatchMaker.cc \
	 LocalMatchMaker.cc  \
	 scoringFunctions.cc \
     MatchedChunk.cc \
	 MapReader.cc \
     xmlWriter.cc \
     FragDatabase.cc \
     seedTypes.cc

_DEPS= AlignmentParams.h\
	array2d.h\
	ContigMapData.h\
	debugUtils.h\
	dp.h\
	exception.h\
	FragDatabase.h\
	GlobalScorer.h\
	globals.h\
	itoa.h\
	localAlignmentScore.h\
	LocalMatchMaker.h\
	LocalScorer.h\
	MapCoord.h\
	MapData.h\
	MapReader.h\
	mapTypes.h\
	MatchedChunk.h\
	match.h\
	MatchMaker.h\
	MatchResult.h\
	MersenneTwister.h\
	OpticalMapData.h\
	parseArgs.h\
	permutation.h\
	ScoreMatrix.h\
	Scorer.h\
	scoringFunctions.h\
	seedTypes.h\
	SeqCoord.h\
	StandardMatchMaker.h\
	types.h\
	utils.h\
	xmlWriter.h\

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJS = $(patsubst %.cc,$(OBJDIR)/%.o,$(srcs))
OBJS_DEBUG = $(patsubst %.cc,$(OBJDIR_DEBUG)/%.o,$(srcs))
OBJS_RELEASE = $(patsubst %.cc,$(OBJDIR_RELEASE)/%.o,$(srcs))


$(OBJS_DEBUG) : $(OBJDIR_DEBUG)/%.o: %.cc | $(OBJDIR_DEBUG)
	@echo making debug object $@
	$(CXX) $(CXXFLAGS) -c $^ -o $@

$(OBJS_RELEASE) : $(OBJDIR_RELEASE)/%.o: %.cc | $(OBJDIR_RELEASE)
	@echo making release object $@
	$(CXX) $(CXXFLAGS) -c $^ -o $@
	

####################################################
# Make archive
../lib/release/soma.a: $(OBJS_RELEASE) | $(LIBDIR_RELEASE)
	@echo making release archive
	ar rvs $@ $(OBJS_RELEASE)

../lib/debug/soma.a: $(OBJS_DEBUG) | $(LIBDIR_DEBUG)
	@echo making debug archive
	ar rvs $@ $(OBJS_DEBUG)

############################################
# Rule to create directories
../build/release ../build/debug: | ../build
	mkdir $@

../lib/release ../lib/debug: | ../lib
	mkdir $@

../bin/release ../bin/debug: | ../bin
	mkdir $@

../test/bin: | ../test
	mkdir $@

../lib ../bin ../build ../test:
	mkdir $@

############################################
# Make executables
.SECONDEXPANSION:
#release debug: main.cc ../lib/$$@/soma.a | ../build/$$@ ../lib/$$@ ../bin/$$@
	#$(CXX) $(CXXFLAGS) main.cc ../lib/$@/soma.a  -o ../bin/$@/match

../bin/release/match: main.cc ../lib/release/soma.a | ../build/release ../lib/release ../bin/release
	$(CXX) $(CXXFLAGS) main.cc ../lib/release/soma.a  -o ../bin/release/match

../bin/debug/match: main.cc ../lib/debug/soma.a | ../build/debug ../lib/debug ../bin/debug
	$(CXX) $(CXXFLAGS) main.cc ../lib/debug/soma.a  -o ../bin/debug/match
   
# TO DO: Fix compilation issues 
../test/bin/testDatabase: ../test/testDatabase.cc ../lib/debug/soma.a | ../test/bin
	$(CXX) $(CXXFLAGS) -I../test -I. ../test/testDatabase.cc $(LIBDIR_DEBUG)/soma.a -o ../test/bin/testDatabase

.PHONY : clean
clean:
	@echo Cleaning Files.
	rm -f ../bin/match ../bin/match_d
	rm -rf ../build
	rm -rf ../lib
