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
#include <core/utils/LinearInterpolator.h>
#include <math.h>
#include <algorithm>

using std::vector;

LinearInterpolator::LinearInterpolator(vector<double> x,
                                       vector<double> y)
{
    n    = (int)x.size();
    mm   = 2;
    jsav = 0;
    cor  = 0;
    xx   = x;
    yy   = y;
    dj = std::max(1, static_cast<int>(std::pow(static_cast<double>(n), 0.25)));
    m_isConfigured = checkConfiguration();
}

LinearInterpolator::~LinearInterpolator()
{

}

double LinearInterpolator::interpolate(double xq)
{
    int jlo = cor ? hunt(xq) : locate(xq);
    return baseInterpolate(jlo, xq);
}

vector<double> LinearInterpolator::interpolateV(vector<double> xq)
{
    vector<double> yq;
    int sz = (int)xq.size();
    yq.resize(sz);
    for(int i=0; i<sz; i++) {
        yq[i] = interpolate(xq[i]);
    }

    return yq;
}

void LinearInterpolator::defineXY(vector<double> x, vector<double> y)
{
    if(x.size() == y.size()) {
        xx = x;
        yy = y;
        n  = (int)x.size();
        jsav = 0;
        cor  = 0;
        dj = std::max(1, static_cast<int>(std::pow(static_cast<double>(n), 0.25)));
        m_isConfigured = true;
    }

}

bool LinearInterpolator::isConfigured()
{
    return m_isConfigured;
}

int LinearInterpolator::locate(const double x)
{
    int ju, jm, jl;
    bool ascnd=(xx[n-1] >= xx[0]);
    jl=0;
    ju=n-1;
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n-mm, jl-((mm-2)>>1)));
}

int LinearInterpolator::hunt(const double x)
{
    int jl=jsav, jm, ju, inc=1;
    bool ascnd=(xx[n-1] >= xx[0]);
    if (jl < 0 || jl > n-1) {
        jl=0;
        ju=n-1;
    } else {
        if (x >= xx[jl] == ascnd) {
            for (;;) {
                ju = jl + inc;
                if (ju >= n-1) { ju = n-1; break;}
                else if (x < xx[ju] == ascnd) break;
                else {
                    jl = ju;
                    inc += inc;
                }
            }
        } else {
            ju = jl;
            for (;;) {
                jl = jl - inc;
                if (jl <= 0) { jl = 0; break;}
                else if (x >= xx[jl] == ascnd) break;
                else {
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n-mm, jl-((mm-2) >> 1)));
}

bool LinearInterpolator::checkConfiguration()
{
    bool status = false;
    if(n < 2) {
        status = false;
    }

    return status;
}

double LinearInterpolator::baseInterpolate(int j, double x)
{
    if (xx[j]==xx[j+1]) {
        return yy[j];
    }

    return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}
