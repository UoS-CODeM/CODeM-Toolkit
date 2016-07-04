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
#ifndef CODEMOPERATORS_H
#define CODEMOPERATORS_H

#include <codemglobal.h>
#include <vector>

namespace CODeM {

// Relations between UncertaintyKernel properties and uncertainty parameters
double linearDecrease(double val);
double skewedIncrease(double val, double alpha);
double skewedDecrease(double val, double alpha);
double lowOnValue(double val, double zeroVal, double width);
double highOnValue(double val, double oneVal, double width);

std::vector<double> directionPerturbation(const std::vector<double> &oVec,
                                          double maxRadius, double pNorm=2);

} // namespace CODeM

#endif // CODEMOPERATORS_H
