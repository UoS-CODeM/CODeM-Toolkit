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
