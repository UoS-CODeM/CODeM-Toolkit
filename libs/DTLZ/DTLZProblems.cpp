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

#include <libs/DTLZ/DTLZProblems.h>
#include <core/CODeMGlobal.h>

#include <math.h>

using std::vector;

namespace DTLZ {

std::vector<double> DTLZ1Modified(const std::vector<double> &x, const int M)
{
    int n = (int)x.size();
    int k = n - M + 1;
    double g = 0.0;
    for (int i = n - k; i < n; i++) {
        g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * CODeM::PI * (x[i] - 0.5));
    }
    // This is the DTLZ paper version, but the huge scaling has no added value
    g = ((double)k + g);

    vector<double> y(M, (1.0 + g) * 0.5);

    for (int i = 0; i < M; i++) {
        int aux = M - (i + 1);
        for (int j = 0; j < aux; j++) {
            y[i] *= x.at(j);
        }
        if (i != 0){
            y[i] *= (1 - x[aux]);
        }
    }

    return y;
}

vector<double > DTLZ1(const vector<double >& x, const int M)
{
    int n = (int)x.size();
    int k = n - M + 1;
    double g = 0.0;
    for (int i = n - k; i < n; i++) {
        g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * CODeM::PI * (x[i] - 0.5));
    }
    g = 100.0 * ((double)k + g);

    vector<double> y(M, (1.0 + g) * 0.5);

    for (int i = 0; i < M; i++) {
        int aux = M - (i + 1);
        for (int j = 0; j < aux; j++) {
            y[i] *= x.at(j);
        }
        if (i != 0){
            y[i] *= (1 - x[aux]);
        }
    }

    return y;
}

vector<double > DTLZ2(const vector<double >& x, const int M)
{
    int i,j;
    int n = (int)x.size();
    int k = n - M + 1;
    double g = 0.0;
    double coss, sine;
    for (int i = n - k; i < n; i++) {
        g += (x[i] - 0.5)*(x[i] - 0.5);
    }

    vector<double> y(M, 0.0);

    for (j=(M-1); j >= 0; j--) {
        coss = 1.0;
        for (i=0; i<M-j-1; i++) {
            coss *= cos(x[i]*CODeM::PI/2.0) ;
        }
        sine = (j>0) ? ((j==M-1) ? sin(x[M-j-1]*CODeM::PI/2.0) : sin(x[i]*CODeM::PI/2.0)) : 1.0;
        y[j] = (1.0+g) * coss * sine;
    }

    return y;
}

}  // namespace DTLZ
