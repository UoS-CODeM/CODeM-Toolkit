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
#include <core/CODeMGlobal.h>

#include <math.h>

using namespace std;

namespace CODeM {

vector<vector<double> > simplexLattice(int h, int k, double s)
{
    vector<vector<double> > wSet;
    //terminating criterion
    if(k==1) {
        wSet.push_back(vector<double>(1, s));
    } else {
        double w = 0.0;
        for(int i = 0; i < h + 1; ++i) {
            vector<vector<double> > kMinusOneSimplex
                    = simplexLattice(h - i, k - 1, s - w);
            for(int j = 0; j < kMinusOneSimplex.size(); ++j) {
                vector<double> ww(kMinusOneSimplex.size() + 1);
                ww[0] = w;
                copy(kMinusOneSimplex[j].begin(), kMinusOneSimplex[j].end(),
                     ww.begin()+1);
                wSet.push_back(ww);
            }
            w += s / h;
        }
    }
    return wSet;
}

} // namespace CODeM
