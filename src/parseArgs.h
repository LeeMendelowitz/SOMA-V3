#include <argp.h>
#include "globals.h"
#include <map>
#include <utility>
#include <string>
#include <iostream>


namespace Parser {

// Function for parsing a single argument
error_t parse_opt(int key, char *arg, struct argp_state *state);

// Options for argp parser
const char * argp_program_version = "SOMA 1.0";
const char * argp_program_bug_address = "<LMendelo@umiacs.umd.edu>";
const char *  doc = "Contig Optical Map Alignment Tool: Align contigs to optical map";
const char * args_doc = "ContigMapFile OpticalMapFile1 [OpticalMapFile 2....]";

// To add a command line argument:
// 1. Add entry to options
// 2. Add entry to ArgParser::printArgs
// 3. Add entry to parse_opt function

// Note: the key values are all initially set to 0, and are subsequently set
typedef struct argp_option argpOption;
argpOption options[] = {
	{"output", 0, "prefix", 0, "prefix for output files. Default: 'output'",0},
    {"circular", 0, 0, 0, "chromosome is circular. Default is non-circular.",0},
//    {"allowFalseCuts", 0, 0, 0, "Allow false cuts in optical map",0},
    {"noReverse", 0, 0, 0, "Do not match reverse of contig to optical map",0},
    {"oneToOneMatch", 0, 0, 0, "Match the entire contig to the entire optical map",0},
//	{"matchFragmentOnce", 0, 0, 0, "greedily match at most one contig to an optical fragment. Default is to allow multiple matches",0},
//	{"allowGaps", 0, 0, 0, "Allow gaps in alignments.",0},
	{"sdMin", 0, "val", 0, "Initial value for the standard deviation threshold. Default: 4",0},
	{"sdMax", 0, "val", 0, "Final value for the standard deviation threshold. Default: 12",0},
	{"pvalue", 0, "threshold", 0, "p-value threshold for permutation test. Default: 0.05",0},
    {"siteCostContig", 0, "value", 0, "cost for a single unsed contig restriction site in an alignment. Default: 36.0",0},
    {"siteCostOptical", 0, "value", 0, "cost for a single unsed optical map restriction site in an alignment. Default: 144.0",0},
    {"maxMissRateContig", 0, "value", 0, "Maximum allowed fraction of contig restriction sites to be unused in an alignment. Default: 0.1",0},
    //{"falseCutRate", 0, "value", 0, "Average number of false cuts per bp in optical map. Note: If provided, allowFalseCuts must be provided.",0},
    {"numPermutationTrials", 0, "value", 0, "Number of trials for permutation test. Default: 0 (i.e. no test)",0},
    {"numThreads", 0, "value", 0, "Number of threads. Default: 1",0},
    //{"maxGapSize", 0, "value", 0, "Maximum open gap size (in bp) Default: 0",0},
	{0}
};

struct argp myArgp = {options, parse_opt, args_doc, doc};


// Wrapper class for the argp parser
class ArgParser
{
    public:
    static ArgParser * instance()
    {
        if (!pInstance_) pInstance_ =  new ArgParser();
        return pInstance_;
    }
    ~ArgParser() { if(pInstance_) delete pInstance_; }
    void parseArgs(int argc, char ** argv);
    void printArgs();
    int getKey(const string& option);
    int getKey(const char * option);
    void getOption(int key, string& option);

    private:
    // private attributes
    static ArgParser * pInstance_;
    map<string, int> optionToKeyMap;
    map<int, string> keyToOptionMap;

    // private methods
    ArgParser(); // Constructor
};

// Create the argp parser
ArgParser::ArgParser()
{
    // Nothing needs to be done. All the action is in parseAgs
}

void ArgParser::parseArgs(int argc, char ** argv)
{
    // Assign the keys for the program options, and populate the optionToKeyMap and keyToOptionMap
    int numOptions = sizeof(options)/sizeof(argp_option);
    for (int i=0; i<numOptions; i++)
    {
        if (!options[i].name) continue;
        options[i].key = i + 1; // Set the key
        string optionName = string(options[i].name);
        optionToKeyMap.insert(make_pair(optionName, options[i].key));
        keyToOptionMap.insert(make_pair(options[i].key, optionName));
    }

    // parse the arguments
    argp_parse(&myArgp, argc, argv, 0, 0, 0);
}

// Print the values of the options
void ArgParser::printArgs()
{

    ostringstream oss;
    int numOpticalMaps = Options::opticalMapList.size();
    for (int i=0; i < numOpticalMaps; i++)
    {
        oss << Options::opticalMapList[i];
        if (i!=(numOpticalMaps-1)) oss << ", ";
    }
    string opMapString = oss.str();

    std::cout << "**********************************************\n"
              << "Running with Inputs:\n\n" 
              << "OpticalMaps: " << opMapString << "\n"
              << "SilicoMap: " << Options::silicoMap << "\n"
              << "Circular: " << Options::circular << "\n"
 //             << "AllowFalseCuts: " << Options::allowFalseCuts << "\n"
//              << "matchFragmentOnce: " << Options::matchFragmentOnce << "\n"
              << "noReverse: " << Options::noReverse << "\n"
              << "oneToOneMatch: " << Options::oneToOneMatch << "\n"
//              << "allowGaps: " << Options::allowGaps << "\n"
              << "sdMin: " << Options::sdMin << "\n"
              << "sdMax: " << Options::sdMax << "\n"
  //            << "maxGapSize: " << Options::maxGapSize << "\n"
              << "pThreshold: " << Options::pThreshold << "\n"
              << "outputPrefix: " << Options::outputPrefix << "\n"
              //<< "C_s: " << Options::C_s << "\n"
              << "C_r_contig: " << Options::C_r_contig << "\n"
              << "C_r_optical: " << Options::C_r_optical << "\n"
   //           << "falseCutRate: " << Options::falseCutRate << "\n"
              << "maxMissRateContig: " << Options::maxMissRateContig << "\n"
              << "numPermutationTrials: " << Options::numPermutationTrials << "\n"
              << "numThreads: " << Options::numThreads << "\n"
              << "**********************************************\n";
}


 /* Parse a single option. */
error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    ArgParser * ap = ArgParser::instance();
    int arg_num = state->arg_num;
    if (key == ap->getKey("circular") ){
        Options::circular = true;
    //} else if (key == ap->getKey("allowFalseCuts") ) {
     //   Options::allowFalseCuts = true;
    //} else if (key == ap->getKey("allowGaps") ) {
        //Options::allowGaps = true;
    //} else if (key == ap->getKey("matchFragmentOnce") ) {
        //Options::matchFragmentOnce = true;
    } else if (key == ap->getKey("noReverse") ) {
        Options::noReverse = true;
    } else if (key == ap->getKey("oneToOneMatch") ) {
        Options::oneToOneMatch = true;
    } else if (key == ap->getKey("sdMin")) {
        Options::sdMin = atoi(arg);
    } else if (key == ap->getKey("sdMax")) {
        Options::sdMax = atoi(arg);
    } else if (key == ap->getKey("pvalue")) {
        Options::pThreshold = atof(arg);
    } else if (key == ap->getKey("output")) {
        Options::outputPrefix = arg;
    } else if (key == ap->getKey("siteCostContig")) {
        Options::C_r_contig = atof(arg);
    } else if (key == ap->getKey("siteCostOptical")) {
        Options::C_r_optical = atof(arg);
    //} else if (key == ap->getKey("falseCutRate")) {
        //Options::falseCutRate = atof(arg);
    } else if (key == ap->getKey("maxMissRateContig")) {
        Options::maxMissRateContig = atof(arg);
    } else if (key == ap->getKey("numPermutationTrials")) {
        Options::numPermutationTrials = atoi(arg);
    } else if (key == ap->getKey("numThreads")) {
        int numThreads = atoi(arg);
        assert(numThreads > 0);
        Options::numThreads = numThreads;
    //} else if (key == ap->getKey("maxGapSize")) {
        //Options::maxGapSize = atoi(arg);
    } else if (key ==  ARGP_KEY_ARG) {
        if (arg_num == 0)
        {
            // Set silico map
            Options::silicoMap = arg;
        }
        else
        {
            // Set optical map
            Options::opticalMapList.push_back(string(arg));
        }
    } else if (key == ARGP_KEY_END) {
        if (arg_num < 2) {
         /* Not enough arguments. */
         argp_usage (state);}
    } else {
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


int ArgParser::getKey(const string& option)
{
    map<string, int>::iterator iter;
    iter = optionToKeyMap.find(option);
    if (iter!= optionToKeyMap.end())
        return iter->second;
    return -1;
}

int ArgParser::getKey(const char * option)
{
    return getKey(string(option));
}

void ArgParser::getOption(int key, string& option)
{
    map<int, string>::iterator iter;
    iter = keyToOptionMap.find(key);
    if (iter != keyToOptionMap.end())
        option =  iter->second;
    option = string("");
}

ArgParser * ArgParser::pInstance_ = 0;

}
