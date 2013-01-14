#include "xmlWriter.h"
#include "MatchResult.h"
#include "utils.h"

#include <cassert>

//
void XMLWriter::open(const char * filename)
{
    // Print the document header
    filename_ = string(filename);
    pXmlFile_ = new ofstream(filename);
    if (pXmlFile_->fail())
    {
        ostringstream msg;
        msg << "ERROR: Could not open file " << filename << " for writing.";
        delete pXmlFile_; pXmlFile_ = 0;
        throw(Exception(msg.str()));
    }
    // Write header
    doc_.save(*pXmlFile_); 
    // Write the root node
    *pXmlFile_ << "<alignments>\n";
}

//
XMLWriter::~XMLWriter()
{
    // close the file
    if(pXmlFile_)
    {
        close();
        delete pXmlFile_;
    }
}

//
void XMLWriter::close()
{
    if (pXmlFile_->good())
        *pXmlFile_ << "</alignments>" << endl;
    pXmlFile_->close();
}

//
void XMLWriter::writeMatchResult(const MatchResult& mr)
{

    xml_node matchResultNode = doc_.append_child("MatchResult");

    // Contig data
    xml_node contigNode = matchResultNode.append_child("contig"); 
    contigNode.append_child("name").append_child(node_pcdata).set_value(mr.contigId_.c_str());
    contigNode.append_child("length").append_child(node_pcdata).set_value(toCStr(mr.contigLength_));

    // Optical Map
    matchResultNode.append_child("chromosome").append_child(node_pcdata).set_value(mr.chromosomeId_.c_str());

    // Alignment statistics
    matchResultNode.append_child("forward").append_child(node_pcdata).set_value(toCStr(mr.forward_));
    matchResultNode.append_child("score").append_child(node_pcdata).set_value(toCStr(mr.score_));
    matchResultNode.append_child("pval").append_child(node_pcdata).set_value(toCStr(mr.pval_));
    matchResultNode.append_child("chi2").append_child(node_pcdata).set_value(toCStr(mr.chi2_));
    matchResultNode.append_child("cStartIndex").append_child(node_pcdata).set_value(toCStr(mr.cStartIndex_));
    matchResultNode.append_child("cEndIndex").append_child(node_pcdata).set_value(toCStr(mr.cEndIndex_));
    matchResultNode.append_child("cStartBp").append_child(node_pcdata).set_value(toCStr(mr.cStartBp_));
    matchResultNode.append_child("cEndBp").append_child(node_pcdata).set_value(toCStr(mr.cEndBp_));
    matchResultNode.append_child("cAlignedBases").append_child(node_pcdata).set_value(toCStr(mr.contigTotalAlignedBases_));
    matchResultNode.append_child("opStartIndex").append_child(node_pcdata).set_value(toCStr(mr.opStartIndex_));
    matchResultNode.append_child("opEndIndex").append_child(node_pcdata).set_value(toCStr(mr.opEndIndex_));
    matchResultNode.append_child("opStartBp").append_child(node_pcdata).set_value(toCStr(mr.opStartBp_));
    matchResultNode.append_child("opEndBp").append_child(node_pcdata).set_value(toCStr(mr.opEndBp_));
    matchResultNode.append_child("opAlignedBases").append_child(node_pcdata).set_value(toCStr(mr.opticalTotalAlignedBases_));
    matchResultNode.append_child("totalHits").append_child(node_pcdata).set_value(toCStr(mr.totalHits_));
    matchResultNode.append_child("totalMisses").append_child(node_pcdata).set_value(toCStr(mr.totalMisses_));
    matchResultNode.append_child("totalMissRate").append_child(node_pcdata).set_value(toCStr(mr.totalMissRate_));
    matchResultNode.append_child("contigHits").append_child(node_pcdata).set_value(toCStr(mr.contigHits_));
    matchResultNode.append_child("contigMisses").append_child(node_pcdata).set_value(toCStr(mr.contigMisses_));
    matchResultNode.append_child("contigMissRate").append_child(node_pcdata).set_value(toCStr(mr.contigMissRate_));
    matchResultNode.append_child("contigUnalignedBases").append_child(node_pcdata).set_value(toCStr(mr.contigUnalignedBases_));
    matchResultNode.append_child("contigUnalignedBaseRatio").append_child(node_pcdata).set_value(toCStr(mr.contigUnalignedBaseRatio_));
    matchResultNode.append_child("contigUnalignedFrags").append_child(node_pcdata).set_value(toCStr(mr.contigUnalignedFrags_));
    matchResultNode.append_child("opticalHits").append_child(node_pcdata).set_value(toCStr(mr.opticalHits_));
    matchResultNode.append_child("opticalMisses").append_child(node_pcdata).set_value(toCStr(mr.opticalMisses_));
    matchResultNode.append_child("opticalMissRate").append_child(node_pcdata).set_value(toCStr(mr.opticalMissRate_));
    matchResultNode.append_child("alignedLengthRatio").append_child(node_pcdata).set_value(toCStr(mr.alignedLengthRatio_));
    matchResultNode.append_child("opticalMatchString").append_child(node_pcdata).set_value(mr.opticalMatchString_.c_str());
    matchResultNode.append_child("opticalAlignedIndex").append_child(node_pcdata).set_value(mr.opticalAlignedIndexString_.c_str());
    matchResultNode.append_child("contigMatchString").append_child(node_pcdata).set_value(mr.contigMatchString_.c_str());
    matchResultNode.append_child("scoreString").append_child(node_pcdata).set_value(mr.scoreString_.c_str());
    matchResultNode.append_child("contigAlignedIndex").append_child(node_pcdata).set_value(mr.contigAlignedIndexString_.c_str());
    matchResultNode.append_child("contigLostIndex").append_child(node_pcdata).set_value(mr.contigLostIndexString_.c_str());


    // Add the matched chunks
    {
    xml_node alignment = matchResultNode.append_child("alignment");
    vector<MatchedChunk>::const_iterator iter = mr.matchedChunkList_.begin();
    const vector<MatchedChunk>::const_iterator E =  mr.matchedChunkList_.end();
    for(; iter != E; iter++)
    {
        addMatchedChunk(alignment, *iter);
    }
    }
    
    // Write the node to the file
    matchResultNode.print(*pXmlFile_);

    // Remove the node from the tree (to save memory)
    doc_.remove_child(matchResultNode);
}

// Write a chunk node
void XMLWriter::addMatchedChunk(xml_node& parent, const MatchedChunk& mc)
{

    xml_node chunkNode = parent.append_child("chunk");

    // Add optical map coords
    xml_node opCoords = chunkNode.append_child("opIndex");
    opCoords.append_child("start").append_child(node_pcdata).set_value(toCStr(mc.getOpticalStartIndex()));
    opCoords.append_child("end").append_child(node_pcdata).set_value(toCStr(mc.getOpticalEndIndex()));

    // Add contig coords
    xml_node contigCoords = chunkNode.append_child("cIndex");
    contigCoords.append_child("start").append_child(node_pcdata).set_value(toCStr(mc.getContigStartIndex()));
    contigCoords.append_child("end").append_child(node_pcdata).set_value(toCStr(mc.getContigEndIndex()));

    // Add optical bp coords
    xml_node opBpCoords = chunkNode.append_child("opBpIndex");
    opBpCoords.append_child("start").append_child(node_pcdata).set_value(toCStr(mc.getOpticalStartBp()));
    opBpCoords.append_child("end").append_child(node_pcdata).set_value(toCStr(mc.getOpticalEndBp()));

    // Add contig bp coords
    xml_node cBpCoords = chunkNode.append_child("cBpCoords");
    cBpCoords.append_child("start").append_child(node_pcdata).set_value(toCStr(mc.getContigStartBp()));
    cBpCoords.append_child("end").append_child(node_pcdata).set_value(toCStr(mc.getContigEndBp()));

    // Add contig frags
    {
    xml_node contigFrags = chunkNode.append_child("contigFrags");
    vector<FragData>::const_iterator iter = mc.getContigFragB();
    vector<FragData>::const_iterator E = mc.getContigFragE();
    for(; iter != E; iter++)
        contigFrags.append_child("frag").append_child(node_pcdata).set_value(toCStr(iter->size_));
    }

    // Add optical frags
    {
    xml_node opticalFrags = chunkNode.append_child("opticalFrags");
    vector<FragData>::const_iterator iter = mc.getOpticalFragB();
    vector<FragData>::const_iterator E = mc.getOpticalFragE();
    for(; iter != E; iter++)
        opticalFrags.append_child("frag").append_child(node_pcdata).set_value(toCStr(iter->size_));
    }

    // Add isContigGap
    chunkNode.append_child("isContigGap").append_child(node_pcdata).set_value(toCStr(mc.isContigGap()));

    // Add isBoundaryChunk
    chunkNode.append_child("isboundarychunk").append_child(node_pcdata).set_value(toCStr(mc.isBoundaryChunk()));

    // Add Score
    const Score& score = mc.getScore();
    xml_node score_node = chunkNode.append_child("score");
    score_node.append_child("contig").append_child(node_pcdata).set_value(toCStr(score.contig));
    score_node.append_child("optical").append_child(node_pcdata).set_value(toCStr(score.optical));
    score_node.append_child("sizing").append_child(node_pcdata).set_value(toCStr(score.sizing));
}
