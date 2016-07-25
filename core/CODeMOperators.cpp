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
#include <core/CODeMOperators.h>
#include <core/utils/ScalingUtils.h>
#include <random>
#include <math.h>

using namespace CODeM::Utils;

namespace CODeM {

double linearDecrease(double val)
{
    return 1.0 - val;
}

double skewedIncrease(double val, double alpha)
{
    return (alpha<0) ? 1.0 : pow(val, alpha);
}

double skewedDecrease(double val, double alpha)
{
    return (alpha<0) ? 0.0 : 1.0 - pow(val, alpha);
}

double lowOnValue(double val, double zeroVal, double width)
{
    double ret = 4.0 / pow(width, 2.0) * pow(val-zeroVal , 2.0);
    return (ret>1.0) ? 1.0 : ret;
}

double highOnValue(double val, double oneVal, double width)
{
    double ret = 1.0 - 4.0 / pow(width, 2.0) * pow(val-oneVal, 2.0);
    return (ret<0.0) ? 0.0 : ret;
}

std::vector<double> directionPerturbation(const std::vector<double> &oVec,
                                     double maxRadius, double pNorm)
{
    // project on the k-1 simplex
    std::vector<double> newObjVec(oVec);
    toUnitVec(newObjVec, 1.0);

    // calculate the p-distance
    double dist = magnitudeP(oVec, pNorm);

    // perturb within a sphere with r=maxRadius
    double s = 0.0;
    for(int i=0; i<newObjVec.size(); i++) {
        double rd = (randUni() * 2.0 -1.0) * sqrt(maxRadius*maxRadius - s);
        s += pow(rd, 2.0);
        newObjVec[i] += rd;
    }

    // project on the p-norm unit sphere
    toUnitVec(newObjVec, pNorm);

    //scale back
    scale(newObjVec, dist);

    return newObjVec;
}

} // namespace CODeM
