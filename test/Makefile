$(CXX)=g++
CXXFLAGS=-Wall

debug: CXX += -DDEBUG -g
debug: debug_executable
SRC=../src

srcs=$(SRC)/match.cc\
     $(SRC)/globals.cc \
     $(SRC)/dp.cc \
     $(SRC)/Scorer.cc \
	 $(SRC)/utils.cc \
	 $(SRC)/debugUtils.cc \
	 $(SRC)/ContigMapData.cc \
	 $(SRC)/MatchResult.cc \
	 $(SRC)/OpticalMapData.cc \
	 $(SRC)/pugixml/pugixml.cpp \
	 $(SRC)/StandardMatchMaker.cc \
	 $(SRC)/LocalMatchMaker.cc  \
	 $(SRC)/scoringFunctions.cc \
     $(SRC)/MatchedChunk.cc \
	 $(SRC)/MapReader.cc \
     $(SRC)/xmlWriter.cc \
	 $(SRC)/FragDatabase.cc \
	 $(SRC)/seedTypes.cc

release: CXX += -O3 -DNDEBUG
release: executable


.cc:	 
	g++ -O0 -g $@.cc -o ../bin/$@

testDatabase: testDatabase.cc
	$(CXX) $(CXXFLAGS) -fopenmp testDatabase.cc $(srcs) -I ../src -o testDatabase

clean:
	rm -f ../bin/match ../bin/schedule *~


