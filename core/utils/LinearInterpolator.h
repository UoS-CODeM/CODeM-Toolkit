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
#ifndef LINEARINTERPOLATOR_H
#define LINEARINTERPOLATOR_H

#include <codemglobal.h>

#include <vector>

namespace CODeM {
namespace Utils {

class LinearInterpolator
{
public:
    LinearInterpolator(std::vector<double> xv, std::vector<double> yv);
    ~LinearInterpolator();

    double interpolate(double xq);
    std::vector<double> interpolateV(std::vector<double> xq);
    virtual void defineXY(std::vector<double> x, std::vector<double> y);
    bool isConfigured();

protected:
    double baseInterpolate(int j, double x);
    int locate(const double x);
    int hunt(const double x);
    virtual bool checkConfiguration();

    int  n;
    int  mm;
    int  jsav;
    int  cor;
    int  dj;
    bool m_isConfigured;
    std::vector<double> xx;
    std::vector<double> yy;
};

} // namespace Utils
} // namespace CODeM
#endif // LINEARINTERPOLATOR_H
