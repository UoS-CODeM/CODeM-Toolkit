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
#ifndef CODEMPROBLEMS_H
#define CODEMPROBLEMS_H

#include <vector>
using std::vector;

namespace CODeM {

vector<double>          CODeM1(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM1(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM1Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM2(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM2(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM2Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM3(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM3(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM3Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM4(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM4(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM4Perturb(const vector<double> &oVec, int nSamp = 1);

// CODeM5Perturb must have both decision and objective vectors defined
vector<double>          CODeM5(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM5(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM5Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM6(const vector<double> &iVec, int nObj);
vector<vector<double> > CODeM6(const vector<double> &iVec, int nObj, int nSamp);
vector<vector<double> > CODeM6Perturb(size_t iVecSize,
                                      const vector<double> &oVec, int nSamp = 1);
vector<vector<double> > CODeM6Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec, int nSamp = 1);

vector<double>          GECCOExample(const vector<double> &iVec, int nObj);
vector<vector<double> > GECCOExample(const vector<double> &iVec, int nObj,
                                     int nSamp);
vector<vector<double> > GECCOExamplePerturb(size_t iVecSize,
                                            const vector<double> &oVec,
                                            int nSamp = 1);
vector<vector<double> > GECCOExamplePerturb(const vector<double> &iVec,
                                            const vector<double> &oVec,
                                            int nSamp = 1);

vector<double> deterministicOVec(int prob, const vector<double> &iVec, int nObj, int k=0);

void createInputBounds(vector<double> &lBounds, vector<double> &uBounds, int prob);

bool validArgs(const vector<vector<double> >          &dVectors,
               const vector<vector<double> >          &oVecDeterm,
               const vector<vector<vector<double> > > &oVecSamps,
               int problem, int k);

void optimalSet(vector<vector<double> >          &dVectors,
                vector<vector<double> >          &oVecDeterm,
                vector<vector<vector<double> > > &oVecSamps,
                int problem, int k);
void randomSet (vector<vector<double> >          &dVectors,
                vector<vector<double> >          &oVecDeterm,
                vector<vector<vector<double> > > &oVecSamps,
                int problem, int k);

} // namespace CODeM

#endif // CODEMPROBLEMS_H
