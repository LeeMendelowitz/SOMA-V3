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
	 MapReader.cc \
     xmlWriter.cc

release: CXX += -O3 -DNDEBUG
release: executable


.cc:	 
	g++ -O0 -g $@.cc -o ../bin/$@

executable: main.cc match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp main.cc $(srcs) -o ../bin/match
    
debug_executable: main.cc match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp main.cc $(srcs) -o ../bin/match_d

testDatabase: test/testDatabase.cc
	$(CXX) $(CXXFLAGS) -fopenmp test/testDatabase.cc $(srcs) -o ../test/testDatabase

clean:
	rm -f ../bin/match ../bin/schedule *~


