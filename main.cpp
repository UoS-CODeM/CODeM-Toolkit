/****************************************************************************
**
** The MIT License (MIT)
**
** Copyright (c) 2016 The University of Sheffield (www.sheffield.ac.uk)
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is
** furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE
**
****************************************************************************/
#include <misc/CODeMMisc.h>
#include <misc/examples/CODeMProblems.h>
#include <core/CODeMGlobal.h>

#include <random>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;
using namespace CODeM;

inline void defineSeed(int seed) {std::srand(seed);}
inline void randomSeed() {std::srand((unsigned)std::time(0));}

void printVector(vector<double> vec, string sep=", ", string endVec="; ")
{
    for(auto i = vec.begin(); i != (vec.end() -1); ++i) {
        cout << *i << sep;
    }
    cout << vec.back() << endVec;
}

void showUsage(char* progName)
{
    cout << "\nUsage: " << progName << " [OPTION(S)]\n\n";

    cout <<
"This program evaluates a set of candidate solutions for an uncertain           \n"
" multiobjective optimization problem using the CODeM toolkit. The output of the\n"
" program is a set of decision vectors and a set of objective vectors for each  \n"
" decision vector. The objective vectors for every solution represent different \n"
" samples from the random variate.                                              \n\n"

"Options:                                                                       \n"
"--------                                                                       \n"
" -h, --help                Print this summary and exit.                        \n\n"
"     --default             Run the default test problem from the GECCO'16 paper\n"
"                           using default settings as specified below.          \n\n"
//" -i, --inFile  = FILENAME  An input file with all the configuration options.   \n"
//"                           If FILENAME is not specified, the options are       \n"
//"                           configured from the command. FILENAME may include a \n"
//"                           relative or absolute path.                          \n\n"
" -f, --file    = FILENAME  A file to write the outputs. If FILENAME is not     \n"
"                           specified, the output is printed to the console.    \n\n"
" -p, --problem = NUMBER    A chioce of benchmark problem from the CODeM suite. \n"
"                           Use NUMBER = 0 for the problem in the GECCO'16      \n"
"                           paper. Use NUMBER = 1,...,6 for CODeM1,...,CODeM6.  \n"
"                           Default is NUMBER = 0.                              \n\n"
//" -x, --solSet  = SET       A set of decision vectors to evaluate. SET needs to \n"
//"                           be provided within double qoutes, where each vector \n"
//"                           is separated with a semicolon, and the elements     \n"
//"                           within a vector separated with a space. e.g., a set \n"
//"                           of three vectors with two variables each:           \n"
//"                           SET = \"1.1 1.2; 2.1 2.2; 3.0 4.0\"                 \n"
//"                           If SET is not specified two defalut sets are        \n"
//"                           generated: one with solutions that are optimal for  \n"
//"                           the deterministic problem, and one with random      \n"
//"                           vectors. Their sizes are specified with the -s and  \n"
//"                           -n options.                                         \n\n"
" -m, --nObj     = NUMBER   The dimensionality of the objective space. If NUMBER\n"
"                           is not specified, the default is 2 objectives.      \n\n"
" -d, --nVars    = NUMBER   The dimensionality of the decision space. Must be   \n"
"                           larger than nObj. If NUMBER is not specified, the   \n"
"                           default is 10.                                      \n\n"
" -s, --nSols    = NUMBER   The number of decision vectors to evaluate. Two     \n"
"                           sets of size NUMBER are generated: one with         \n"
"                           solutions that are optimal for the deterministic    \n"
"                           problem, and one with random vectors. If NUMBER is  \n"
"                           not specified, the default value is NUMBER = 10.    \n\n"
" -n, --nSamps   = NUMBER   The number of function evaluations for each decision\n"
"                           vector. if NUMER is not specified, the default value\n"
"                           is NUMBER = 5.                                      \n\n"
" -r, --rndSeed  = NUMBER   The seed for the pseudo-random number generator. If \n"
"                           NUMBER is not specified, the default seed is 0. For \n"
"                           a random seed, based on the CPU time, provide a     \n"
"                           negative value for NUMBER.                          \n"
" -k, --nDirVars = NUMBER   The number of direction related variables for WFG   \n"
"                           based problems. If NUMBER is not specified, the     \n"
"                           default of nObj - 1 is used                         \n\n";
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        showUsage(argv[0]);
        return EXIT_SUCCESS;
    }

    /// Set default values
    int seed   = 0;
    int nObj   = 2;
    int nVars  = 10;
    int nSols  = 10;
    int nSamps = 5;
    int prob   = 0;
    int k      = 0;


    /// Parse command line inputs
    int argInd = 1;
    while(argInd < argc) {
        string arg = argv[argInd++];

        if ((arg == "-h") || (arg == "--help")) {
            showUsage(argv[0]);
            return EXIT_SUCCESS;

        } else if (arg == "--default") {
            if(argc != 2) {
                cerr << "--default option must be the only argument." << endl;
                return EXIT_FAILURE;
            }// else run with default settings

        } else if ((arg == "-f") || (arg == "--file")) {
            if (argInd < argc) {
                freopen(argv[argInd++], "w", stdout);
            } else {
                cerr << "--file option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-p") || (arg == "--problem")) {
            if (argInd < argc) {
                int argI = atoi(argv[argInd++]);
                if((argI >= 0) && (argI <= 6)) {
                    prob = argI;
                } else {
                    cerr << "Invalid argument for --problem option: Requires a "
                            "number between 0-6." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                cerr << "--problem option requires one argument." << endl;
                return EXIT_FAILURE;
            }
            // TODO: accept a set to evaluate
            //        } else if ((arg == "-x") || (arg == "--solSet")) {
            //            if (argInd < argc) {

            //            } else {
            //                cerr << "--solSet option requires one argument." << endl;
            //                return EXIT_FAILURE;
            //            }

        } else if ((arg == "-m") || (arg == "--nObj")) {
            if (argInd < argc) {
                int argI = atoi(argv[argInd++]);
                if(argI >= 2) {
                    nObj = argI;
                } else {
                    cerr << "Invalid argument for --nObj option: Number of "
                            "objectives must be larger than 1." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                cerr << "--nObj option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-d") || (arg == "--nVars")) {
            if (argInd < argc) {
                int argI = atoi(argv[argInd++]);
                if(argI >= 3) {
                    nVars = argI;
                } else {
                    cerr << "Invalid argument for --nVars option: Number of "
                            "variables must be larger than 2." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                cerr << "--nVars option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-s") || (arg == "--nSols")) {
            if (argInd < argc) {
                int argI = atoi(argv[argInd++]);
                if(argI >= 0) {
                    nSols = argI;
                } else {
                    cerr << "Invalid argument for --nSols option: Number of "
                            "solutions must be larger than 0." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                cerr << "--nSols option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-n") || (arg == "--nSamps")) {
            if (argInd < argc) {
                int argI = atoi(argv[argInd++]);
                if(argI >= 0) {
                    nSamps = argI;
                } else {
                    cerr << "Invalid argument for --nSamps option: Number of "
                            "solutions must be larger than 0." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                cerr << "--nSamps option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-r") || (arg == "--rndSeed")) {
            if (argInd < argc) {
                seed = atoi(argv[argInd++]);
            } else {
                cerr << "--rndSeed option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else if ((arg == "-k") || (arg == "--nDirVars")) {
            if (argInd < argc) {
                k = atoi(argv[argInd++]);
            } else {
                cerr << "--nDirVars option requires one argument." << endl;
                return EXIT_FAILURE;
            }

        } else {
            cerr << "Unknown argument " << arg << endl;
            return EXIT_FAILURE;
        }
    }

    /// Define the random seed
    if(seed < 0) {
        randomSeed();
    } else {
        defineSeed(seed);
    }

    if(k == 0) {
        k = nObj - 1;
    }


    /// Print run configuration settings
    cout << "% CODeM Toolkit Demosntrator v1.0\n"
            "% Copyright (c) 2016 The University of Sheffield\n\n"
         << "% Configuration settings:"    << endl
         << "rndSeed  = " << seed   << ";" << endl
         << "nObj     = " << nObj   << ";" << endl
         << "nVars    = " << nVars  << ";" << endl
         << "nSols    = " << nSols  << ";" << endl
         << "nSamps   = " << nSamps << ";" << endl
         << "problem  = " << prob   << ";" << endl
         << "nDirVars = " << k      << ";" << endl;


    /// Construct the set of solutions
    vector<vector<double> > dVectors(nSols, vector<double>(nVars));
    vector<vector<double> > oVecDeterm(nSols, vector<double>(nObj));
    vector<vector<vector<double> > > oVecSamps(nSols,
                       vector<vector<double> >(nSamps, vector<double>(nObj)));

    /// Assign Pareto optimal values
    if(prob == 5) {
        cout << "\n% Optimal set cannot be analytically derived for CODeM5" << endl;
    }
    else
    {
        optimalSet(dVectors, oVecDeterm, oVecSamps, prob, k);

        // Display the results
        cout << "\n% Optimal decision vectors:" << endl;
        cout << "optSol = [";
        for(int i = 0; i < dVectors.size(); ++i) {
            printVector(dVectors[i]);
        }
        cout << "];" << endl;

        cout << "\n% Optimal deterministic objective vectors:" << endl;
        cout << "optDetermObj = [";
        for(int i = 0; i < oVecDeterm.size(); ++i) {
            printVector(oVecDeterm[i]);
        }
        cout << "];" << endl;

        cout << "\n% Samples for optimal solutions:" << endl;
        for(int v = 0; v < oVecSamps.size(); ++v) {
            cout << "optObjSamps{" << v+1 << "} = [";
            for(int i = 0; i < oVecSamps[v].size(); ++i) {
                printVector(oVecSamps[v][i]);
            }
            cout << "];" << endl;
        }
    }

    /// Random solutions
    randomSet(dVectors, oVecDeterm, oVecSamps, prob, k);

    // Display the results
    cout << "\n% Random decision vectors:" << endl;
    cout << "rndSol = [";
    for(int i = 0; i < dVectors.size(); ++i) {
        printVector(dVectors[i]);
    }
    cout << "];" << endl;

    cout << "\n% Random deterministic objective vectors:" << endl;
    cout << "rndDetermObj = [";
    for(int i = 0; i < oVecDeterm.size(); ++i) {
        printVector(oVecDeterm[i]);
    }
    cout << "];" << endl;

    cout << "\n% Samples for random solutions:" << endl;
    for(int v = 0; v < oVecSamps.size(); ++v) {
        cout << "rndObjSamps{" << v+1 << "} = [";
        for(int i = 0; i < oVecSamps[v].size(); ++i) {
            printVector(oVecSamps[v][i]);
        }
        cout << "];" << endl;
    }

    return EXIT_SUCCESS;
}
