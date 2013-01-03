$(CXX)=g++
CXXFLAGS=-Wall

debug: CXX += -DDEBUG -g
debug: debug_executable

release: CXX += -O3 -DNDEBUG
release: executable


install:
	make globals; make utils; make match; make schedule
.cc:	 
	g++ -O0 -g $@.cc -o ../bin/$@

executable: match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp globals.cc dp.cc match.cc utils.cc MatchResult.cc ContigMapData.cc OpticalMapData.cc pugixml/pugixml.cpp -o ../bin/match
    
debug_executable: match.cc globals.cc dp.cc utils.cc ContigMapData.cc MatchResult.cc OpticalMapData.cc
	$(CXX) $(CXXFLAGS) -fopenmp globals.cc dp.cc match.cc utils.cc MatchResult.cc \
	          ContigMapData.cc \
			  OpticalMapData.cc \
			  pugixml/pugixml.cpp -o ../bin/match_d

clean:
	rm -f ../bin/match ../bin/schedule *~
