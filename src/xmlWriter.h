#ifndef XMLWRITER_H
#define XMLWRITER_H

#include "pugixml.h"
#include <fstream>

using namespace pugi;
using namespace std;


class MatchResult;
class MatchedChunk;

class XMLWriter
{
    public:

        ~XMLWriter();
        void open(const char * filename);
        void close();
        void writeMatchResult(const MatchResult& mr);

    private:
        string filename_;
        xml_document doc_;
        ofstream * pXmlFile_;

        void addMatchedChunk(xml_node& parent, const MatchedChunk& mc);
};


#endif
