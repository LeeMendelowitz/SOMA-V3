$(CXX)=g++
CXXFLAGS=-Wall

debug: CXX += -DDEBUG -g
debug: debug_executable

srcs=match.cc\
     globals.cc \
     dp.cc \
     Scorer.cc \
	 utils.cc \
	 debugUtils.cc \
	 ContigMapData.cc \
	 MatchResult.cc \
	 OpticalMapData.cc \
	 pugixml/pugixml.cpp \
	 StandardMatchMaker.cc \
	 LocalMatchMaker.cc  \
	 scoringFunctions.cc \
     MatchedChunk.cc \


release: CXX += -O3 -DNDEBUG
release: executable


install:
	make globals; make utils; make match; make schedule
.cc:	 
	g++ -O0 -g $@.cc -o ../bin/$@

executable: match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp $(srcs) -o ../bin/match
    
debug_executable: match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp $(srcs) -o ../bin/match_d

clean:
	rm -f ../bin/match ../bin/schedule *~
