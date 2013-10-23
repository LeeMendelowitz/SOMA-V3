CXXFLAGS=-Wall -fopenmp -pedantic -std=c++0x                  
IDIR=

GTEST_INC = /home/lee/packages/gtest-1.7.0/include

OBJDIR_RELEASE = ../build/release
OBJDIR_DEBUG = ../build/debug
OBJDIR_GPROF = ../build/gprof

LIBDIR_RELEASE = ../lib/release
LIBDIR_DEBUG = ../lib/debug
LIBDIR_GPROF = ../lib/gprof

all: debug release gprof testDatabase

debug: CXXFLAGS += -DDEBUG -g -O0
debug: OBJDIR = $(OBJDIR_DEBUG)
debug: LIBDIR = $(LIBDIR_DEBUG)
debug: DIR=debug
debug: ../bin/debug/match

release: CXXFLAGS += -O3 -DNDEBUG
release: OBJDIR = $(OBJDIR_RELEASE)
release: LIBDIR = $(LIBDIR_RELEASE)
release: DIR=release
release: ../bin/release/match

gprof: CXXFLAGS += -DDEBUG -pg -O0
gprof: OBJDIR = $(OBJDIR_DEBUG)
gprof: LIBDIR = $(LIBDIR_DEBUG)
gprof: DIR=gprof
gprof: ../bin/gprof/match



testDatabase: ../test/bin/testDatabase
testDatabase: CXXFLAGS += -DDEBUG -g -O0

gtestDatabase: ../test/bin/gtestDatabase
#gtestDatabase: CXXFLAGS += -DDEBUG -g -O0
gtestDummy: ../test/bin/gtestDummy

srcs=match.cc\
     globals.cc \
     dp.cc \
	 seededDp.cc \
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
	 MapData.cc \
     xmlWriter.cc \
     ChunkDatabase.cc \
     MapChunk.cc \
	 MapChunkUtils.cc

_DEPS= AlignmentParams.h\
	array2d.h\
	ContigMapData.h\
	debugUtils.h\
	dp.h\
	seededDp.h\
	exception.h\
	ChunkDatabase.h\
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
	ScoreMatrixSeeded.h\
	ScoreCell.h\
	ScorePathStep.h\
	Scorer.h\
	scoringFunctions.h\
	MapChunk.h\
	MapChunkUtils.h\
	SeqCoord.h\
	StandardMatchMaker.h\
	types.h\
	utils.h\
	xmlWriter.h\

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJS = $(patsubst %.cc,$(OBJDIR)/%.o,$(srcs))
OBJS_DEBUG = $(patsubst %.cc,$(OBJDIR_DEBUG)/%.o,$(srcs))
OBJS_RELEASE = $(patsubst %.cc,$(OBJDIR_RELEASE)/%.o,$(srcs))
OBJS_GPROF = $(patsubst %.cc,$(OBJDIR_GPROF)/%.o,$(srcs))


%.o : %.cc $(DEPS)

$(OBJS_DEBUG) : $(OBJDIR_DEBUG)/%.o: %.cc | $(OBJDIR_DEBUG)
	@echo making debug object $@
	$(CXX) $(CXXFLAGS) -c $^ -o $@

$(OBJS_RELEASE) : $(OBJDIR_RELEASE)/%.o: %.cc | $(OBJDIR_RELEASE)
	@echo making release object $@
	$(CXX) $(CXXFLAGS) -c $^ -o $@

$(OBJS_GPROF) : $(OBJDIR_GPROF)/%.o: %.cc | $(OBJDIR_GPROF)
	@echo making gprof object $@
	$(CXX) $(CXXFLAGS) -c $^ -o $@
	

####################################################
# Make archive
../lib/release/soma.a: $(OBJS_RELEASE) | $(LIBDIR_RELEASE)
	@echo making release archive
	ar rvs $@ $(OBJS_RELEASE)

../lib/gprof/soma.a: $(OBJS_GPROF) | $(LIBDIR_GPROF)
	@echo making gprof archive
	ar rvs $@ $(OBJS_GPROF)

../lib/debug/soma.a: $(OBJS_DEBUG) | $(LIBDIR_DEBUG)
	@echo making debug archive
	ar rvs $@ $(OBJS_DEBUG)

############################################
# Rule to create directories
../build/release ../build/debug ../build/gprof: | ../build
	mkdir $@

../lib/release ../lib/debug ../lib/gprof: | ../lib
	mkdir $@

../bin/release ../bin/debug ../bin/gprof: | ../bin
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
	@echo "Building release version of SOMA"
	$(CXX) $(CXXFLAGS) main.cc ../lib/release/soma.a  -o ../bin/release/match

../bin/debug/match: main.cc ../lib/debug/soma.a | ../build/debug ../lib/debug ../bin/debug
	$(CXX) $(CXXFLAGS) main.cc ../lib/debug/soma.a  -o ../bin/debug/match

../bin/gprof/match: main.cc ../lib/gprof/soma.a | ../build/gprof ../lib/gprof ../bin/gprof
	$(CXX) $(CXXFLAGS) main.cc ../lib/gprof/soma.a  -o ../bin/gprof/match

   
../test/bin/testDatabase: ../test/testDatabase.cc ../lib/debug/soma.a | ../test/bin
	$(CXX) $(CXXFLAGS) -I../test -I. ../test/testDatabase.cc $(LIBDIR_DEBUG)/soma.a -o $@

#../test/bin/gtestDatabase: ../test/gtest_database.cc ../lib/debug/soma.a | ../test/bin
	#$(CXX) $(CXXFLAGS) -I../test -I. -I$(GTEST_INC) /usr/local/lib/gtest.a ../lib/debug/soma.a -o $@

../test/bin/%: ../test/%.cc ../lib/debug/soma.a | ../test/bin
	$(CXX) $(CXXFLAGS) -lpthread -I../test -I. -I$(GTEST_INC) /usr/local/lib/gtest.a ../lib/debug/soma.a -o $@

.PHONY : clean
clean:
	@echo Cleaning Files.
	rm -f ../bin/match ../bin/match_d
	rm -rf ../build
	rm -rf ../lib
