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
#ifndef ABSTRACTINTERPOLATOR_H
#define ABSTRACTINTERPOLATOR_H




#include <QtMath>
#include <vector>

class AbstractInterpolator
{
public:
    explicit AbstractInterpolator(vector<double> x, vector<double> y, int m);
    virtual ~AbstractInterpolator();

    double interpolate(double xq);
    vector<double> interpolateV(vector<double> xq);
    virtual void defineXY(vector<double> x, vector<double> y);
    bool isConfigured();

protected:
    int locate(const double x);
    int hunt(const double x);
    virtual double baseInterpolate(int jlo, double x) = 0;
    virtual bool checkConfiguration();

    int n;
    int mm;
    int jsav;
    int cor;
    int dj;
    bool m_isConfigured;
    vector<double> xx;
    vector<double> yy;
};

#endif // ABSTRACTINTERPOLATOR_H
