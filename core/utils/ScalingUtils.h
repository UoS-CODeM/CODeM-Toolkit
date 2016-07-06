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
