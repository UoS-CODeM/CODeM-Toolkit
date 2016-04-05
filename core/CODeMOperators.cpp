/****************************************************************************
**
** Copyright (C) 2012-2015 The University of Sheffield (www.sheffield.ac.uk)
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
#include <core/CODeMOperators.h>
#include <random>
#include <tigon/Utils/NormalisationUtils.h>
#include <qmath.h>

namespace CODeM {

double linearDecrease(double val)
{
    return 1.0 - val;
}

double skewedIncrease(double val, double alpha)
{
    return (alpha<0) ? 1.0 : qPow(val, alpha);
}

double skewedDecrease(double val, double alpha)
{
    return (alpha<0) ? 0.0 : 1.0 - qPow(val, alpha);
}

double lowOnValue(double val, double zeroVal, double width)
{
    double ret = 4.0 / qPow(width, 2.0) * qPow(val-zeroVal , 2.0);
    return (ret>1.0) ? 1.0 : ret;
}

double highOnValue(double val, double oneVal, double width)
{
    double ret = 1.0 - 4.0 / qPow(width, 2.0) * qPow(val-oneVal, 2.0);
    return (ret<0.0) ? 0.0 : ret;
}

vector<double> directionPerturbation(const vector<double> oVec,
                                     double maxRadius, double pNorm)
{
    // project on the k-1 simplex
    vector<double> newObjVec(oVec);
    toUnitVec(newObjVec, 1.0);

    // calculate the p-distance
    vector<double> vec(oVec);
    double dist = magnitudeAndDirectionP(vec, pNorm);

    // perturb within a sphere with r=maxRadius
    double s = 0.0;
    for(int i=0; i<newObjVec.size(); i++) {
        double rd = TRAND.randUni(2.0, -1.0) *
                qSqrt(maxRadius*maxRadius - s);
        s += qPow(rd, 2.0);
        newObjVec[i] += rd;
    }

    // project on the p-norm unit sphere
    toUnitVec(newObjVec, pNorm);

    //scale back
    scale(newObjVec, dist);

    return newObjVec;
}

} // namespace CODeM
