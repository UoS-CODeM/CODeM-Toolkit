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
#include <core/utils/ScalingUtils.h>

#include <math.h>

using namespace std;

namespace CODeM {
namespace Utils {

double magnitudeP(const vector<double> &vec, double p)
{
    double magnitude = 0.0;
    for(int i=0; i < vec.size(); i++) {
        magnitude += pow(vec[i], p);
    }
    magnitude = pow(magnitude, 1.0 / p);

    return magnitude;
}

void toUnitVec(vector<double>& vec, double norm)
{
    double magnitude = magnitudeP(vec, norm);
    scale(vec, 1.0/magnitude);
}

double distanceP(const vector<double> &m, const vector<double> &c, double p)
{
    double dist = 0.0;
    if(m.size() == c.size()) {
        for(int i=0; i<m.size(); i++) {
            dist += pow(abs(m[i] - c[i]), p);
        }
        dist = pow(dist, 1.0 / p);
        return dist;
    } else {
        return -1.0;
    }
}

void scale(vector<double>& vec, double factor)
{
    for(int i=0; i < vec.size(); i++) {
        vec[i] *= factor;
    }
}

void normaliseToUnitInterval(double& val, double  lBound, double  uBound)
{
    if(val <= lBound) {
        val = 0.0;
    } else if(val >= uBound) {
            val = 1.0;
    } else {
        val = (val - lBound) / (uBound - lBound);
    }
}

void scaleBackFromUnitInterval(double& normVal, double  lBound, double  uBound)
{
    if(normVal >= 1) {
        normVal = uBound;
    } else if(normVal <= 0.0) {
        normVal = lBound;
    } else {
        normVal = lBound + normVal * (uBound - lBound);
    }
}

void normaliseToUnitBox(vector<double>& vec,
                        const vector<double> &lBounds, const vector<double> &uBounds)
{
    for(int i=0; i < vec.size(); i++) {
        if(vec[i] <= lBounds[i]) {
            vec[i] = 0.0;
        }else if(vec[i] >= uBounds[i]) {
            vec[i] = 1.0;
        }   else {
            vec[i] = (vec[i] - lBounds[i]) / (uBounds[i] - lBounds[i]);
        }
    }
}


void scaleBackFromUnitBox(vector<double>& vec,
                          const vector<double> &lBounds, const vector<double> &uBounds)
{
    for(int i=0; i < vec.size(); i++) {
        if(vec[i] > 1.0) {
            vec[i] = uBounds[i] - (1.0 - vec[i]) * (uBounds[i] - lBounds[i]);
        } else if(vec[i] < 0.0) {
            vec[i] = lBounds[i] - vec[i] * (uBounds[i] - lBounds[i]);
        } else {
            vec[i] = lBounds[i] + vec[i] * (uBounds[i] - lBounds[i]);
        }
    }
}

double directedBoxedIntervalLength(const vector<double> &dir, double p)
{
    double maxComponent = 0.0;
    for(int i=0; i < dir.size(); i++) {
        if(dir[i] > maxComponent) {
            maxComponent = dir[i];
        }
    }
    double len = 0.0;
    for(int i=0; i < dir.size(); i++) {
        len += pow((dir[i] / maxComponent), p);
    }
    len = pow(len, 1.0 / p);
    return len;
}

} // namespace Utils
} // namespace CODeM
