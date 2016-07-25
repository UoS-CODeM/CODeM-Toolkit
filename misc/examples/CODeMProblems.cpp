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
#include <misc/examples/CODeMProblems.h>
#include <core/UncertaintyKernel.h>
#include <core/CODeMOperators.h>
#include <core/CODeMDistribution.h>
#include <core/RandomDistributions.h>
#include <core/utils/ScalingUtils.h>

#include <libs/WFG/ExampleProblems.h>
#include <libs/DTLZ/DTLZProblems.h>

using std::vector;
using namespace WFGT::Toolkit::Examples::Problems;
using namespace CODeM::Utils;

namespace CODeM {

vector<double> CODeM1(const vector<double> &iVec, int k, int nObj)

{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM1Perturb(oVec)[0];
}

vector<vector<double> > CODeM1(const vector<double> &iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM1Perturb(oVec, nSamp);
}

vector<vector<double> > CODeM1Perturb(const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0);
    vector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double peakTend, peakLoc, dirPertRad;

    peakTend = uk.proximity();
    peakLoc  = lowOnValue(uk.proximity(), 0.0, 0.05);

    double distanceNorm = 2.0;

    dirPertRad = 0.0;

    // Create the CODeM distribution
    PeakDistribution* d = new PeakDistribution(peakTend, peakLoc);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> CODeM2(const vector<double> &iVec, int k, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM2Perturb(oVec)[0];
}

vector<vector<double> > CODeM2(const vector<double> &iVec, int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM2Perturb(oVec, nSamp);
}

vector<vector<double> > CODeM2Perturb(const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0.0);
    vector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, dirPertRad;

    uniLB  = uk.proximity();
    uniUB  = uniLB;

    double distanceNorm = 2.0;

    dirPertRad = 0.1 * uk.symmetry();

    // Create the CODeM distribution
    UniformDistribution* d =new UniformDistribution(uniLB, uniUB);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> CODeM3(const vector<double> &iVec, int k, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM3Perturb(oVec)[0];
}

vector<vector<double> > CODeM3(const vector<double> &iVec, int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG4(iVec, k, nObj);

    return CODeM3Perturb(oVec, nSamp);
}

vector<vector<double> > CODeM3Perturb(const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0);
    vector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;


    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniLoc, uniUB, peakTend, peakLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = skewedDecrease(uk.proximity(), 1.5);
    uniUB  = uniLB + (1-uniLoc) * (1.0-uniLB);

    peakTend = uk.proximity();
    peakLoc  = uk.symmetry();

    double distanceNorm = 2.0;

    dirPertRad = 0.04*lowOnValue(uk.oComponent(0), 0.45, 0.3);

    // Create the CODeM distribution
    MergedDistribution* d = new MergedDistribution;
    d->appendDistribution(new UniformDistribution(uniLB, uniUB), 0.5);
    d->appendDistribution(new PeakDistribution(peakTend, peakLoc), 0.5);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> CODeM4(const vector<double> &iVec, int k, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG6(iVec, k, nObj);

    return CODeM4Perturb(oVec)[0];
}

vector<vector<double> > CODeM4(const vector<double> &iVec, int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG6(iVec, k, nObj);

    return CODeM4Perturb(oVec, nSamp);
}

vector<vector<double> > CODeM4Perturb(const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0.0);
    vector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double peakTend, peakLoc, dirPertRad;

    peakTend = uk.proximity()+0.1;
    peakLoc  = 0.8;

    double distanceNorm = 2.0;

    dirPertRad = 0.2*linearDecrease(uk.symmetry())+0.01;

    // Create the CODeM distribution
    PeakDistribution* d = new PeakDistribution(peakTend, peakLoc);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> CODeM5(const vector<double> &iVec, int k, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG8(iVec, k, nObj);

    return CODeM5Perturb(iVec, oVec)[0];
}

vector<vector<double> > CODeM5(const vector<double> &iVec, int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = WFG8(iVec, k, nObj);

    return CODeM5Perturb(iVec, oVec, nSamp);
}

vector<vector<double> > CODeM5Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0.0);
    vector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 4.0*(i+1);
    }
    double lb = 2.0/4.0;
    double ub = 1.0;

    vector<double> inLB(iVec.size());
    vector<double> inUB(iVec.size());

    createInputBounds(inLB, inUB, 5);

    UncertaintyKernel uk(iVec, oVec, lb, ub, inLB, inUB, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, uniLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = linearDecrease(uk.dComponent(0));
    uniUB  = uniLB + (1.0-uniLoc) * (1.0-uniLB);

    double distanceNorm = 2.0;

    dirPertRad = 0.1 * uk.dComponent(0);

    // Create the CODeM distribution
    UniformDistribution* d =new UniformDistribution(uniLB, uniUB);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> CODeM6(const vector<double> &iVec, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = DTLZ::DTLZ1Modified(iVec, nObj);

    return CODeM6Perturb(iVec.size(), oVec)[0];
}

vector<vector<double> > CODeM6(const vector<double> &iVec,
                                int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = DTLZ::DTLZ1Modified(iVec, nObj);

    return CODeM6Perturb(iVec.size(), oVec, nSamp);
}

vector<vector<double> > CODeM6Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec,
                                       int nSamp)
{
    return CODeM6Perturb(iVec.size(), oVec, nSamp);
}

vector<vector<double> > CODeM6Perturb(size_t iVecSize, const vector<double> &oVec, int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0.0);

    // DTLZ1 is modified so the 100 scale of the distance function
    // is not included
    double maxVal = 1.125 * iVecSize;
    vector<double> antiIdeal(oVec.size(), maxVal);

    vector<double> normVec(oVec);
    toUnitVec(normVec, 2.0);
    double sFactor = magnitudeP(normVec, 1);

    double ub = 1.0 / sFactor;
    double lb = 0.5 / maxVal / sFactor;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, uniLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = linearDecrease(uk.proximity()*uk.symmetry());
    uniUB  = uniLB + (1.0-uniLoc) * (1.0-uniLB);

    double distanceNorm = 1.0;

    dirPertRad = 0.2 * uk.oComponent(0);

    // Create the CODeM distribution
    UniformDistribution* d = new UniformDistribution(uniLB, uniUB);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<double> GECCOExample(const vector<double> &iVec, int nObj)
{
    // Evaluate the decision vector
    vector<double> oVec = DTLZ::DTLZ1Modified(iVec, nObj);

    return GECCOExamplePerturb(iVec.size(), oVec)[0];
}

vector<vector<double> > GECCOExample(const vector<double> &iVec, int nObj, int nSamp)
{
    // Evaluate the decision vector
    vector<double> oVec = DTLZ::DTLZ1Modified(iVec, nObj);

    return GECCOExamplePerturb(iVec.size(), oVec, nSamp);
}

vector<vector<double> > GECCOExamplePerturb(size_t iVecSize, const vector<double> &oVec,
                                            int nSamp)
{
    // Set the uncertainty kernel
    vector<double> ideal(oVec.size(), 0.0);

    double maxVal = 1.125 * (iVecSize - oVec.size() + 1) + 0.5;
    vector<double> antiIdeal(oVec.size(), maxVal);

    vector<double> normVec(oVec);
    toUnitVec(normVec, 2.0);
    double sFactor = magnitudeP(normVec, 1);

    double ub = 1.0 / sFactor;
    double lb = 0.5 / maxVal / sFactor;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, uniLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = 0.9 + 0.1 * lowOnValue(uk.proximity(), 0.0, 1.0);
    uniUB  = 1.0 -  uniLoc * (1.0-uniLB);

    double distanceNorm = 1.0;

    dirPertRad = 0.02 + 0.1 * linearDecrease(uk.symmetry());

    // Create the CODeM distribution
    UniformDistribution* d = new UniformDistribution(uniLB, uniUB);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    vector<vector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.push_back(cd.sampleDistribution());
    }
    return samples;
}

vector<vector<double> > GECCOExamplePerturb(const vector<double> &iVec,
                                            const vector<double> &oVec, int nSamp)
{
    return GECCOExamplePerturb(iVec.size(), oVec, nSamp);
}


vector<double> deterministicOVec(int prob, const vector<double> &iVec, int nObj, int k)
{
    vector<double> oVec;

    switch(prob) {
    case 1: case 2: case 3:
        oVec = WFG4(iVec, k, nObj);
        break;

    case 4:
        oVec = WFG6(iVec, k, nObj);
        break;

    case 5:
        oVec = WFG8(iVec, k, nObj);
        break;

    case 6: case 0:
        oVec = DTLZ::DTLZ1Modified(iVec, nObj);
        break;
    default:
        break;
    }

    return oVec;
}

void createInputBounds(vector<double> &lBounds,
                       vector<double> &uBounds,
                       int prob)
{
    size_t nVar = lBounds.size();

    switch(prob) {

    case 1: case 2: case 3: case 4: case 5:
        for(int i = 0; i < nVar; i++) {
            lBounds[i] = 0.0;
            uBounds[i] = 2.0 * (i + 1.0);
        }
        break;

    case 6: case 7:
    default:
        for(int i = 0; i < nVar; i++) {
            lBounds[i] = 0.0;
            uBounds[i] = 1.0;
        }
        break;
    }
}

bool validArgs(const vector<vector<double> >          &dVectors,
               const vector<vector<double> >          &oVecDeterm,
               const vector<vector<vector<double> > > &oVecSamps,
               int problem, int k)
{
    if((problem < 0) || (problem > 6)) {
        return false;
    }

    int nSols = dVectors.size();
    if((nSols < 1) || (oVecDeterm.size() != nSols)
                   || (oVecSamps.size()  != nSols)) {
        return false;
    }

    int nSamps  = oVecSamps[0].size();
    if(nSamps < 1) {
        return false;
    }

    int nVars = dVectors[0].size();
    int nObj  = oVecDeterm[0].size();
    if(oVecSamps[0][0].size() != nObj);

    if((problem == 0) || (problem == 6)) {
        return (nVars > nObj);
    } else {
        return (k >= 1) && (k < nVars) && (nObj >= 2) && (k % ( nObj-1 ) == 0);
    }

}

void optimalSet(vector<vector<double> >          &dVectors,
                vector<vector<double> >          &oVecDeterm,
                vector<vector<vector<double> > > &oVecSamps,
                int problem, int k)
{
    if((problem == 5) ||
            !validArgs(dVectors, oVecDeterm, oVecSamps, problem, k)) {
        return;
    }

    int nSols   = dVectors.size();
    int nVars   = dVectors[0].size();
    int nObj    = oVecDeterm[0].size();
    int nSamps  = oVecSamps[0].size();

    // create optimal decision vectors with random direction variables
    switch(problem)
    {
    case 1: case 2: case 3: case 4:
    {
        for(int i = 0; i < nSols; ++i) {
            for(int j = 0; j < k; ++j) {
                dVectors[i][j] = randUni() * 2.0 * (j + 1);
            }
            for(int j = k; j < nVars; ++j) {
                dVectors[i][j] = 0.35 * 2.0 * (j + 1);
            }
        }
        break;
    }
    case 0: case 6: default:
    {
        k = nObj - 1;
        for(int i = 0; i < nSols; ++i) {
            for(int j = 0; j < k; ++j) {
                dVectors[i][j] = randUni();
            }
            for(int j = k; j < nVars; ++j) {
                dVectors[i][j] = 0.5;
            }
        }
        break;
    }
    }

    // Evaluate the vectors
    for(int i = 0; i < nSols; ++i) {
        oVecDeterm[i] = deterministicOVec(problem, dVectors[i], nObj, k);
        switch(problem)
        {
        case 0: default:
        {
            oVecSamps[i] = GECCOExamplePerturb(nVars, oVecDeterm[i], nSamps);
            break;
        }
        case 1:
        {
            oVecSamps[i] = CODeM1Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 2:
        {
            oVecSamps[i] =  CODeM2Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 3:
        {
            oVecSamps[i] =  CODeM3Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 4:
        {
            oVecSamps[i] =  CODeM4Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 6:
        {
            oVecSamps[i] =  CODeM6Perturb(nVars, oVecDeterm[i], nSamps);
            break;
        }
        }
    }
}

void randomSet(vector<vector<double> >          &dVectors,
               vector<vector<double> >          &oVecDeterm,
               vector<vector<vector<double> > > &oVecSamps,
               int problem, int k)
{
    if(!validArgs(dVectors, oVecDeterm, oVecSamps, problem, k)) {
        return;
    }

    int nSols   = dVectors.size();
    int nVars   = dVectors[0].size();
    int nObj    = oVecDeterm[0].size();
    int nSamps  = oVecSamps[0].size();

    // create random decision vectors
    switch(problem)
    {
    case 1: case 2: case 3: case 4:
    {
        for(int i = 0; i < nSols; ++i) {
            for(int j = 0; j < nVars; ++j) {
                dVectors[i][j] = randUni() * 2.0 * (j + 1);
            }
        }
        break;
    }
    case 0: case 6: default:
    {
        k = nObj - 1;
        for(int i = 0; i < nSols; ++i) {
            for(int j = 0; j < nVars; ++j) {
                dVectors[i][j] = randUni();
            }
        }
        break;
    }
    }

    // Evaluate the vectors
    for(int i = 0; i < nSols; ++i) {
        oVecDeterm[i] = deterministicOVec(problem, dVectors[i], nObj, k);
        switch(problem)
        {
        case 0: default:
        {
            oVecSamps[i] = GECCOExamplePerturb(nVars, oVecDeterm[i], nSamps);
            break;
        }
        case 1:
        {
            oVecSamps[i] = CODeM1Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 2:
        {
            oVecSamps[i] =  CODeM2Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 3:
        {
            oVecSamps[i] =  CODeM3Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 4:
        {
            oVecSamps[i] =  CODeM4Perturb(oVecDeterm[i], nSamps);
            break;
        }
        case 6:
        {
            oVecSamps[i] =  CODeM6Perturb(nVars, oVecDeterm[i], nSamps);
            break;
        }
        }
    }
}

} // namespace CODeM
