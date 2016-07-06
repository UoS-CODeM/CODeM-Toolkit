/****************************************************************************
**
** Copyright (C) 2012-2016 The University of Sheffield (www.sheffield.ac.uk)
**
** This file is part of Liger.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General
** Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
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
    magnitudeP(vec, norm);
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
