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
#ifndef CODEMPROBLEMS_H
#define CODEMPROBLEMS_H

#include <codemglobal.h>

#include <vector>
using std::vector;

namespace CODeM {

vector<double>          CODeM1(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM1(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM1Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM2(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM2(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM2Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM3(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM3(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM3Perturb(const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM4(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM4(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM4Perturb(const vector<double> &oVec, int nSamp = 1);

// CODeM5Perturb must have both decision and objective vectors defined
vector<double>          CODeM5(const vector<double> &iVec, int k, int nObj);
vector<vector<double> > CODeM5(const vector<double> &iVec, int k, int nObj, int nSamp);
vector<vector<double> > CODeM5Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec, int nSamp = 1);

vector<double>          CODeM6(const vector<double> &iVec, int nObj);
vector<vector<double> > CODeM6(const vector<double> &iVec, int nObj, int nSamp);
vector<vector<double> > CODeM6Perturb(const vector<double> &iVec,
                                      const vector<double> &oVec, int nSamp = 1);

vector<double> deterministicOVec(int prob, const vector<double> &iVec, int nObj, int k=0);

void createInputBounds(vector<double> &lBounds, vector<double> &uBounds, int prob);
} // namespace CODeM

#endif // CODEMPROBLEMS_H
