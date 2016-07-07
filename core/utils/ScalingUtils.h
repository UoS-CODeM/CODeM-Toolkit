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
#ifndef NORMALISATIONUTILS_H
#define NORMALISATIONUTILS_H

#include <codemglobal.h>
#include <vector>

using std::vector;
namespace CODeM {
namespace Utils {

// Returns the p-norm of the vector
double magnitudeP(const vector<double>& vec, double p = 2.0);

void toUnitVec(vector<double>& vec, double norm = 2.0);

// The p-norm of the difference between two vectors
double distanceP(const vector<double> &m, const vector<double> &c, double p=2.0);

// The vector is scaled by a factor
void scale(vector<double>& dir, double factor);

// The vector is normlaised to be within a unit hyperbox,
// where the ideal and the anti-ideal are set to 0 and 1, respectively.
void normaliseToUnitInterval(double& val, double  lBound, double  uBound);

void scaleBackFromUnitInterval(double& normVal, double  lBound, double  uBound);

void normaliseToUnitBox(vector<double>& vec,
                        const vector<double> &lBounds, const vector<double> &uBounds);

void scaleBackFromUnitBox(vector<double>& normVec,
                          const vector<double> &lBounds, const vector<double> &uBounds);

// The length of a vector passing through the 1-norm direction vector dir,
// ending at the boundary of the unit box.
double directedBoxedIntervalLength(const vector<double> &dir, double p=2.0);

} // namespace Utils
} // namespace CODeM
#endif // NORMALISATIONUTILS_H
