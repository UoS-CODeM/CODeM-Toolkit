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
#include <core/utils/LinearInterpolator.h>
#include <math.h>
#include <algorithm>

using std::vector;

namespace CODeM {
namespace Utils {

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

} // namespace Utils
} // namespace CODeM
