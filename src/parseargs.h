#ifndef PARSEARGS_H
#define PARSEARGS_H
#include <iostream>
#include <unistd.h>
#include "globals.h"
#include <argp.h>

using namespace std;
using namespace Globals; // globals in namespace are set by ParseArgs

void PrintUsage(const char* s) {
    cerr << "\nUSAGE: " << s << "  [-s in-silico-map] [-m optical-map] [-c circular-flag] [-r match once]"
         << "[-p permutation-test-threshold] [-f F-test-threshold] [-o output-prefix]\n\n";
    return;
}

void PrintHelp(const char* s) {
    PrintUsage (s);
    cerr
            << "-s path          Set path of in-silico map file to use\n"
            << "-m path          Set path of optical map file to use\n"
            << "-c flag          1 if map is circular and 0 otherwise\n"
            << "-r               Match a fragment from optical map to at most one contig\n"
            << "-p threshold     P-value threshold for the permutation test\n"
            << "-f threshold     P-value threshold for the F-test\n"
            << "-o path          Prefix for output files\n"
            << "-h               Display help information\n"
            << endl;
    cerr
            << "Aligns in-silico map with optical map\n"
            << endl;
    return;
}

// ParseArgs will read the commmand line arguments and set global variables corresponding to
// the inputs
void ParseArgs(int argc, char** argv) {
    int ch, errflg = 0;
    optarg = NULL;
    while(!errflg  && ((ch = getopt(argc, argv, "s:m:c:r:p:f:o:h")) != EOF))
    {
        switch (ch) {
        case 's':
            OPT_SilicoMap = optarg; // path of in-silico map
            break;
        case 'm':
            OPT_OpticalMap = optarg; // path of optical map
            break;
        case 'c':
            OPT_CircularFlag = atoi(optarg); // 1 if map is circular
            break;
        case 'r':
            OPT_MatchOnce = true;
            break;
        case 'p':
            OPT_PT_threshold = atof(optarg); // P-value threshold for permutation test
            break;
        case 'f':
            OPT_FT_threshold = atof(optarg); // P-value threshold for the F-test
            break;
        case 'o':
            OPT_OutputPrefix = optarg; // prefix for output files
            break;
        case 'h':
            PrintHelp(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        default:
            errflg++;
        }
    }

    if(argc == 1 || errflg > 0 || optind != argc) {
        PrintUsage(argv[0]);
        cerr << "Try '" << argv[0] << " -h' for more information.\n";
        exit(EXIT_FAILURE);
    }
}
#endif
