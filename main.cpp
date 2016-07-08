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
#include <string>

using namespace std;
using namespace CODeM;

inline void defineSeed(int seed) {std::srand(seed);}
inline void randomSeed() {std::srand((unsigned)std::time(0));}

void dispVector(vector<double> vec, string sep=", ", string endVec="; ")
{
    for(auto i = vec.begin(); i != (vec.end() -1); ++i) {
        cout << *i << sep;
    }
    cout << vec.back() << endVec;
}

int main()
{
    defineSeed(0);

    int nObj = 2;
    int nVar = 6;
    int nPareto = 5;
    int nRand   = 10;
    int nSamp = 50;

    vector<vector<double> > paretoDeterministic;
    vector<vector<double> > randDeterministic;
    vector<vector<vector<double> > > paretoObj;
    vector<vector<vector<double> > > randObj;

    vector<vector<double> > paretoDirVars = simplexLattice(nPareto-1, 2);

    // Pareto optimal vectors
    for(int i=0; i<paretoDirVars.size(); i++) {
        vector<double> iVec(nVar, 0.5);
        iVec[0] = paretoDirVars[i][0];
        vector<double> determObjVec = deterministicOVec(7, iVec, nObj);
        paretoDeterministic.push_back(determObjVec);
        vector<vector<double> > oVecSamps = CODeM::GECCOExample(iVec, nObj, nSamp);
        paretoObj.push_back(oVecSamps);
    }

    // Random vectors
    for(int i=0; i<nRand; i++) {
        vector<double> iVec;
        for(int j=0; j<nVar; j++) {
            iVec.push_back(randUni());
        }
        vector<double> determObjVec = deterministicOVec(7, iVec, nObj);
        randDeterministic.push_back(determObjVec);
        vector<vector<double> > oVecSamps = CODeM::GECCOExample(iVec, nObj, nSamp);
        randObj.push_back(oVecSamps);
    }

    // Display the results
    cout << "\n% Optimal vectors:" << endl;
    for(int v=0; v<paretoObj.size(); v++) {
        cout << "paretoObj{" << v+1 << "} = [";
        for(int i=0; i<paretoObj[v].size(); i++) {
            dispVector(paretoObj[v][i]);
        }
        cout << "];" << endl;
    }

    cout << "\n% Random vectors:" << endl;
    for(int v=0; v<randObj.size(); v++) {
        cout << "randObj{" << v+1 << "} = [";
        for(int i=0; i<randObj[v].size(); i++) {
            dispVector(randObj[v][i]);
        }
        cout << "];" << endl;
    }

    cout << "\n% Optimal deterministic vectors:" << endl;
    cout << "determOptimal = [";
    for(int i=0; i<paretoDeterministic.size(); i++) {
        dispVector(paretoDeterministic[i]);
    }
    cout << "];" << endl;

    cout << "\n% Random deterministic vectors:" << endl;
    cout << "determRand = [";
    for(int i=0; i<randDeterministic.size(); i++) {
        dispVector(randDeterministic[i]);
    }
    cout << "];\n" << endl;

    return 0;
}
